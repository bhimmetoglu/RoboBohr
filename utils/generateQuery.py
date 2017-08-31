# Author    : Burak Himmetoglu
# Date in   : 04-13-2017
# Date curr : 08-10-2017
# -- AlcNet -- # 

import numpy as np
import os
import re
import pickle
import gzip
import json
import io
import pandas as pd
from tqdm import tqdm
from collections import Counter
from scipy.spatial.distance import pdist, squareform
from Constants import *

class molecule(object):
  """ A class for storing molecular data """
  def __init__ (self, sid):
    self.sid = sid # Id number of the molecule

  def atoms(self, positions, names):
    """ Atoms in the molecule """
    self.positions = positions # Atomic positions
    self.names = names         # Names

  def properties(self,mmff94_en, shape_multipoles, self_overlap):
    """ Features extracted from sdf files """
    self.mmff94_en = mmff94_en
    self.shape_multipoles = shape_multipoles
    self.self_overlap = self_overlap 

  def returnDict(self):
    """ Returns a dictionary for storing in json """
    return_dict = {}

    return_dict['id'] = self.sid                           # Id
    return_dict['En'] = self.mmff94_en                     # Energy
    return_dict['shapeM'] = self.shape_multipoles.tolist() # Shape multipoles
    return_dict['self_overlap'] = self.self_overlap        # Self overlap 
      
    list_atoms = []                                         # Atoms
    for el in zip(self.names, self.positions):
      d = {}
      d['type'] = el[0]
      d['xyz'] = el[1].tolist()
      list_atoms.append(d)

    return_dict['atoms'] = list_atoms

    # Return 
    return return_dict

class pubChemQuery(object):
  def __init__ (self, sdfFile):
    self.sdfFile = sdfFile # Path to the sdf file

  def readData(self):
    # Check that the SDF folder exists 
    assert os.path.exists(self.sdfFile) == True, "Error: SDF file does not exist"
        
    # Read the sdf files in the list and concatenate to fParse
    f = gzip.open(self.sdfFile,"r")
    fParse = f.readlines()
    f.close()

    # Decode
    fParse = [item.decode("utf-8") for item in fParse]    

    # Save in object
    self.fParse = fParse

  def count(self):
    """ Count the atoms and find the maximum """
    natMax = 0 # initiate
    for i, line in enumerate(self.fParse):
      if re.search("-OEChem-", line):
        nat = int(self.fParse[i+2][0:3])
        natMax = max(natMax, nat)

    # Return maximum number of atoms
    return natMax

  def processQuery(self, periodic_table, return_dictionary = True, natMax = 50):
    """ Process the query and save input files """ 
    auToA = 0.529 # au is 0.529 A

    # Read molecule sdf's
    startMolecule = []  # Line where molecule data starts
    ids = []            # Id of the molecule
    mmff94 = []         # MMFF94 energies
    shp_multp = []      # Shape multipoles 
    f_so = []           # Self overlaps (?)

    # Go over all the lines
    print("Reading sdf file...\n")
    ind = 0
    for line in tqdm(self.fParse):
      if re.search("-OEChem-", line):
        startMolecule.append(ind+2)
      if re.search("<PUBCHEM_COMPOUND_CID>", line):
        ids.append(ind+1)
      if re.search("<PUBCHEM_MMFF94_ENERGY>", line):
        mmff94.append(ind+1)
      if re.search("<PUBCHEM_SHAPE_MULTIPOLES>", line):
        shp_multp.append(ind+1)
      if re.search("<PUBCHEM_FEATURE_SELFOVERLAP>", line):
        f_so.append(ind+1)

      ind +=1

    print("Done!\n")

    # For each molecule read the structure
    print("Constructing molecule objects...\n")
    listOfMolecules = []
    ind = 0
    elementSet = set(periodic_table.keys())
    for imol in tqdm(startMolecule):
      # Initiate a temporary molecule
      ID = int(self.fParse[ids[ind]].split()[0])
      tempMolecule = molecule(sid = ID)    

      # Initiate
      index = startMolecule[ind]
      natoms = int(self.fParse[index][0:3])
      endCoords = index + natoms

      # Initiate lists for storing molecule info
      x = []; y = []; z = []; elements = []
      for iat in range(index+1,endCoords+1):
        text = self.fParse[iat].split()
        x.append(text[0]); y.append(text[1]); z.append(text[2])
        elements.append(text[3])

      # Store features
      mmff94_en = float(self.fParse[mmff94[ind]].split()[0])
      self_overlap = float(self.fParse[f_so[ind]].split()[0])
      shape_multipoles = np.zeros(14)
      for m in range(14):
        shape_multipoles[m] = float(self.fParse[shp_multp[ind]+m].split()[0])

      tempMolecule.properties(mmff94_en, shape_multipoles, self_overlap)

      # Find unique terms in element and count 
      typ = list(Counter(elements).keys())
      ntype = len(typ)

      # If there are elements that are not included in the periodic table, 
      # or the molecule contains more than natMax atoms skip 
      if ( not (elementSet > set(typ)) or natoms > natMax):
        ind += 1
        continue

      # Save positions in a np.array then compute maximum distance between the atoms
      positions = np.column_stack( (np.array(x, dtype=float),
                                    np.array(y, dtype=float),
                                    np.array(z, dtype=float)) )

      # Save in tempMolecule
      tempMolecule.atoms(positions = positions, names = elements)

      # Append to List of Molecules
      if (return_dictionary):
        molecule_data = tempMolecule.returnDict()
        listOfMolecules.append(molecule_data)
      else:
        listOfMolecules.append(tempMolecule)

      # Iterate
      ind += 1

    # Print report
    print ("Constructed {:} molecules".format(len(listOfMolecules)))

    # return listOfMolecules
    return listOfMolecules


# Function for storing JSON data
def storeJSON(listOfMolecules, filename):
  """ Store JSON """
  print("Saving {:d} molecules...\n".format(len(listOfMolecules)))

  try:
    to_unicode = unicode
  except NameError:
    to_unicode = str
    
  with io.open(filename, 'w', encoding='utf-8') as outfile:
    json_data = json.dumps(listOfMolecules, indent=4, sort_keys=True, 
                           separators=(',', ': '), ensure_ascii=False)
    outfile.write(to_unicode(json_data))

  print("Done!\n")


#### ---- Generation of Coulomnb Matrices ---- ####

# Function to shift to center of mass and rescale
def rescale_molecule(molecule, periodic_table):
  """ Shift to cm and rescale """
  elements = molecule.names
  pos = np.array(molecule.positions, dtype = float)
  
  # Store masses
  natoms = len(elements)
  masses = np.zeros(natoms)
  for at in range(natoms):
    masses[at] = periodic_table[molecule.names[at]][1]

  # Center of mass coordinates
  posCm = np.dot(masses, pos) / np.sum(masses)
    
  # Rescale (to au)
  posCm  = (pos - posCm)/auToA

  return posCm

# Function that converts listOfMolecules into input data
def storeMolecules(listOfMolecules, natMax, start_index = 0, end_index = 1):
  """ Function for writing data on disk """
  
  # Number of molecules in the list
  numMolecules = len(listOfMolecules)

  # Initialize a list
  print("Saving processed molecular data...\n")

  # Initiate
  data_CM = np.zeros((numMolecules, natMax*(natMax+1)//2), dtype=float)
  data_ids = np.zeros(numMolecules, dtype=int)
  data_multipoles = np.zeros((numMolecules,14), dtype=float)
  data_overlap = np.zeros(numMolecules, dtype=float)
  data_mmff94 = np.zeros(numMolecules, dtype=float)

  for ind in tqdm(range(numMolecules)):
    tempMolecule = listOfMolecules[ind]
    natoms = len(tempMolecule.names)   # Number of atoms in the molecule
    pubchem_id = int(tempMolecule.sid) # PubChem Id
    pos = rescale_molecule(tempMolecule, periodic_table) # Rescaled positions
 
    # Construct Coulomb Matrices
    tiny = 1e-20
    dm = pdist(pos)

    # Z's
    Z = np.array([periodic_table[tempMolecule.names[at]][0] for at in range(natoms)],dtype=float)

    # Coulomb matrix 
    coulomb_matrix = np.outer(Z,Z) / (squareform(dm) + tiny)

    # Full CM padded with zeroes
    full_CM = np.zeros((natMax, natMax))
    full_Z = np.zeros(natMax)
    full_CM[0:natoms, 0:natoms] = coulomb_matrix 
    full_Z[0:natoms] = Z

    # Coulomb vector
    iu = np.triu_indices(natMax,k=1)  # No diagonal k=1
    coulomb_vector = full_CM[iu]

    # Order
    shuffle = np.argsort(-coulomb_vector)
    coulomb_vector = coulomb_vector[shuffle]

    # Construct feature vector
    coulomb_matrix = squareform(coulomb_vector)
    assert np.trace(coulomb_matrix) == 0, "Wrong Coulomb Matrix!"

    coulomb_matrix += 0.5*np.power(full_Z,2.4)*np.eye(natMax) # Add the diagonal terms
    iu = np.triu_indices(natMax) # Upper diagonal 
    feature_vector = coulomb_matrix[iu]
    assert feature_vector.shape[0] == natMax*(natMax+1)//2, "Wrong feature dimensions"

    # Save data 
    data_CM[ind] = feature_vector
    data_multipoles[ind] = tempMolecule.shape_multipoles
    data_overlap[ind] = tempMolecule.self_overlap
    data_mmff94[ind] = tempMolecule.mmff94_en
    data_ids[ind] = pubchem_id

    # Index for saving
    if (ind == 0):
       start_index = pubchem_id
    elif (ind == (numMolecules-1)):
       end_index = pubchem_id

  # Pandas data frame
  dat = np.column_stack((data_CM, data_multipoles, data_overlap))
  df = pd.DataFrame(dat)

  # Column names
  numfeats = np.shape(dat)[1]
  cols = [x for x in range(1,numfeats+1,1)]
  col_names = list(map(lambda x: 'f'+str(x), cols))
  df.columns = col_names

  # Add mmff94 and id
  df.insert(0, 'pubchem_id', data_ids)
  df['En'] = data_mmff94

  # Save
  path = "molecules_" + str(start_index) + "_" + str(end_index) + ".csv"
  df.to_csv(path)

  print("Done!\n")

# Coulomb vector
#iu = np.triu_indices(natMax,k=1)  # No diagonal k=1
#coulomb_vector = full_CM[iu]

# Order
#shuffle = np.argsort(-coulomb_vector)
#coulomb_vector = coulomb_vector[shuffle]
