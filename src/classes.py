# Author: Burak Himmetoglu
# Date  : 08-18-2016
# -- Project RoboBohr -- # 

import numpy as np
from collections import Counter
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from Constants import *
import pairFeatures 

class atom:
  """ Atom class: Not currently used """
  def __init__ (self, position, typ, name, mass):
    self.position  # Coordinates x,y,z
    self.typ = typ # Type of the atom
    self.name = name  # Name of the atom (i.e. the name of the pseudo-potential)
    self.mass = mass # Atomic mass

class molecule:
  """ Molecule class """
  def __init__ (self, sid):
    self.sid = sid # Id number of the molecule

  def atoms(self, positions, names, elInd):
    """ Atoms in the molecule """
    self.positions = positions # Atomic positions
    self.names = names # Names
    self.elInd = elInd # Index corressponding to the element

class featureMatrix:
  """ The Coulomb Matrix """
  # Initiate with an instance of molecule
  def __init__ (self, molecule):
    self.molecule = molecule

  def indexAtoms(self,masses):
    """ Shift atoms to center of mass coordinates, then reorder; 
        Atoms closest to center of mass has index 0, etc
    """
    # Get the molecular structure
    elements = self.molecule.names # Names of elements
    pos = np.array(self.molecule.positions, dtype = float) # Atomic positions  
    typ = Counter(elements).keys() # Atomic types, i.e. unique names
    ntype = len(typ) # Number of distinct types  
    natoms = len(elements) # Total number of atoms in the molecule

    # Match masses
    m = []
    for iat in range(natoms):
      elementIndex = elementList.index(elements[iat])
      m.append(masses[elementIndex])

    # Center of mass coordinates
    m = np.array(m, dtype = float)
    rCm = np.dot(m, pos) / np.sum(m)

    # Shift coordinates by rCm
    pos += -rCm

    # Distances to cm
    distCm = np.apply_along_axis(np.linalg.norm,1,pos)
    increasingOrder = np.argsort(distCm)

    # Reorder the molecule with increasing order and cm shifted positions
    elements = np.array(elements)
    self.molecule.names = elements[increasingOrder]
    self.molecule.positions = pos[increasingOrder,:]

  def coulombMatrix(self,elementList, ZList, Zpow = 2.4, eigenval = True, nrandom = 0):
    """ Construction of Coulomb Matrices (or eigenspectrum) """

    # Get the molecular structure 
    elements = self.molecule.names # Names of elements
    pos = np.array(self.molecule.positions, dtype = float) # Atomic positions  
    typ = Counter(elements).keys() # Atomic types, i.e. unique names
    ntype = len(typ) # Number of distinct types  
    natoms = len(elements) # Total number of atoms in the molecule

    # Match Z's
    Z = []
    for iat in range(natoms):
      elementIndex = elementList.index(elements[iat])
      Z.append(ZList[elementIndex])

    # Convert positions to au
    auToA = 0.529 # au is 0.529 A
    pos /= auToA

    # Convert Z's into an array
    Z = np.array(Z, dtype = float)

    # Compute the Coulomb Matrix
    CM = 1/(squareform(pdist(pos,"euclidean")) + 1e-10); np.fill_diagonal(CM,0)
    CM = np.outer(Z,Z) * CM

    # Add diagonal terms
    CM += 0.5 * Z**(Zpow) * np.eye(natoms)

    # Return the Coulomb Matrix
    if (eigenval == True):
      # Now, we want to compute the eigenvalues, and use them as features
      w,v = np.linalg.eig(CM)
      w = np.sort(w)
      return w

    elif (eigenval == False and nrandom == 0):
      # Return the full Coulomb Matrix
      return CM

    elif (eigenval == False and nrandom != 0):
      # Initiate a new matrix for storage
      CMr = np.zeros((natoms,natoms,nrandom))
      # Compute norm of rows
      Cn = np.apply_along_axis(np.linalg.norm, 1, CM)
      for ir in range(nrandom):
        # Random eps
        eps = np.random.randn(natoms)
        Ctemp = Cn + eps
        # Pertumation
        P = Ctemp.argsort()
        CMr[:,:,ir] = (CM[:,P])[P,:]
      # Return CMr
      return CMr

  def pairFeatureMatrix(self, elementList):
    """ Construction of pair-distance matrices """

    # Initiate
    nSpecies = len(elementList)
    
    # Get the molecular structure 
    pos = np.array(self.molecule.positions, dtype = float) # Atomic positions  
    elInd = np.array(self.molecule.elInd, dtype = np.intc) # Element indices matching to elementList
    natoms = len(self.molecule.names) # Total number of atoms in the molecule
 
    # Initiate the matrix
    dim1 = natoms * (natoms -1)/2 # First dimension (pairwise distances)
    dim2 = nSpecies * (nSpecies + 1)/2 # Number of possible pairs
    featMat = np.zeros((dim1,dim2)) # To be passed to fun_pairFeatures (compiled C code)

    # Call the C function to store the pairFeatures
    pairFeatures.fun_pairFeatures(nSpecies, natoms, elInd, pos, featMat)

    # Return featMat
    return featMat
