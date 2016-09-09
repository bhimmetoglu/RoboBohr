# Author: Burak Himmetoglu
# Date  : 08-18-2016
# -- Project RoboBohr -- # 

import numpy as np
import os
import re
from collections import Counter
from scipy.spatial.distance import pdist
from classes import *
from Constants import *

class pubChemQuery:
  def __init__ (self, sdfFolder):
    self.sdfFolder = sdfFolder # Path to the folder where SDF files reside

  def readData(self):
    # Initial check
    if (not os.path.exists(self.sdfFolder)):
      print("Error: SDF folder does not exist!")
      quit()
        
    # List of the sdf files in df folder
    sdfList = os.listdir(self.sdfFolder)

    # Read the sdf files in the list and concatenate to fParse
    fParse = []
    for filSdf in sdfList:
      filnam = os.path.join(self.sdfFolder,filSdf)
      f = open(filnam,"r")
      fParse += f.readlines()
      f.close()

    # Save in object
    self.fParse = fParse

  def processQuery(self, elementList, valences, natMax, box, tol = 5.29):
    """ Process the query and save input files 
        elementList : List of elements we want to include
    """
    auToA = 0.529 # au is 0.529 A
    # Read molecule sdf's
    startMolecule = []; ids = []
    ind = 0
    for line in self.fParse:
      if re.search("-OEChem-", line):
        startMolecule.append(ind+2)
      if re.search("<PUBCHEM_COMPOUND_CID>", line):
        ids.append(ind+1)

      ind +=1

    # For each molecule read the structure
    listOfMolecules = []
    ind = 0
    for imol in startMolecule:
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

      # Find unique terms in element and count 
      typ = Counter(elements).keys()
      typCount = Counter(elements).values(); ntype = len(typ)

      # If the elements are in the list of elements_list, total number of atoms <= nat_limit, 
      # keep going; Otherwise proceed with next index
      if ( not (set(elementList) > set(typ)) or natoms > natMax or natoms < 2):
        ind += 1
        continue

      # Make sure that the molecule is closed shell, pass otherwise
      nelec = 0
      for it in range(ntype):
        elementIndex = elementList.index(typ[it])
        nelec += typCount[it] * valences[elementIndex]
      if ( nelec % 2 != 0):
        ind += 1
        continue

      # Save positions in a np.array then compute maximum distance between the atoms
      positions = np.column_stack( (np.array(x),np.array(y),np.array(z)) )

      #print index, typ
      distMatrix = pdist(positions,"euclidean"); maxDist = np.max(distMatrix)
      # If the molecule fits into the box, keep going
      if ( box * auToA - maxDist  <= tol):
         ind +=1
         continue

      # Create the element Index (elInd)
      elInd = []
      for iat in range(natoms):
        elementIndex = elementList.index(elements[iat])
        elInd.append(elementIndex)

      # Save in tempMolecule
      tempMolecule.atoms(positions = positions, names = elements, elInd = np.array(elInd))

      # Append to List of Molecules
      listOfMolecules.append(tempMolecule)

      # Iterate
      ind += 1

    # return listOfMolecules
    return listOfMolecules
