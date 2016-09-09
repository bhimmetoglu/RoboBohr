# Author: Burak Himmetoglu
# Date  : 08-18-2016
# -- Project RoboBohr -- # 

import os
import numpy as np
from glob import glob
from collections import Counter
from Constants import *
from classes import *

def storeFeatures(listOfMolecules, natMax = 50, eigenval = True, nrandom = 0):
  """ A function for writing Coulomb Matrices in file """

  # Initiate data matrices
  numMolecules = len(listOfMolecules)
  if (eigenval == True):
    Xdata = np.zeros((numMolecules,natMax+1))
  elif (eigenval == False and nrandom == 0):
    nTot = natMax * (natMax + 1)/2; iu = np.triu_indices(natMax) # Upper triangular part
    Xdata = np.zeros((numMolecules,nTot+1))
  else:
    nTot = natMax * (natMax + 1)/2; iu = np.triu_indices(natMax) # Upper triangular part
    Xdata = np.zeros((numMolecules*nrandom,nTot+1))

  # Loop over listOfMolecules
  for ind in range(numMolecules):
    tempMolecule = listOfMolecules[ind]
    natoms = len(tempMolecule.names)
    # Create an instance of featureMatrix
    fM = featureMatrix(molecule = tempMolecule)
    # Reorder atoms with respect to center of mass coordinates (if desired)
    # fM.indexAtoms(masses)
    # Get the Coulomb Matrix
    CM = fM.coulombMatrix(elementList, ZList, Zpow = 2.4, eigenval = eigenval, nrandom = nrandom)

    if (eigenval == True):
      Xdata[ind,0] = tempMolecule.sid
      Xdata[ind,1:natoms+1] = CM
    if (eigenval == False and nrandom == 0):
      XTemp = np.zeros((natMax,natMax))
      XTemp[:natoms,:natoms] = CM
      Xdata[ind,0] = tempMolecule.sid
      Xdata[ind,1:] = XTemp[iu] #Upper triangular part
    if (eigenval == False and nrandom !=0):
      XTemp = np.zeros((natMax,natMax))
      for ir in range(nrandom):
        ind_ir = ir + nrandom * ind
        XTemp[:natoms,:natoms] = CM[:,:,ir]
        Xdata[ind_ir,0] = tempMolecule.sid
        Xdata[ind_ir,1:] = XTemp[iu] # Upper triangular part

  # Write features in file
  if (not os.path.exists(pathData)):
    os.mkdir(pathData)
  
  filNam = "coulombX" + ".csv"
  filX = os.path.join(pathData,filNam)
  np.savetxt(filX,Xdata,delimiter=",",fmt="%1.6f")

def storePairFeatures(listOfMolecules, natMax = 50):
  """ Uses the distance matrices as features """

  # Initiate data matrices
  nSpecies = len(elementList)
  numMolecules = len(listOfMolecules)
  dim1 = (natMax * (natMax-1)/2); dim2 = (nSpecies * (nSpecies+1)/2) 
  dimFeat = dim1 * dim2

  # +1 for id number
  Xdata = np.zeros((numMolecules,dimFeat+1))

  # Loop over listOfMolecules
  for ind in range(numMolecules):
    tempMolecule = listOfMolecules[ind]
    natoms = len(tempMolecule.names)
    # Create an instance of featureMatrix
    fM = featureMatrix(molecule = tempMolecule)
    # Get the pair-distance matrix
    pDM = fM.pairFeatureMatrix(elementList)

    # Add pDM to Xdata
    XTemp = np.zeros((dim1,dim2))
    tempDim1 = (natoms * (natoms-1)/2)
    XTemp[:tempDim1,:] = pDM
    Xdata[ind,0] = tempMolecule.sid
    Xdata[ind,1:] = XTemp.reshape(dimFeat) # Unroll 

  # Write features in file
  if (not os.path.exists(pathData)):
    os.mkdir(pathData)

  filNam = "pairFeatMat" + ".csv"
  filX = os.path.join(pathData,filNam)
  np.savetxt(filX,Xdata,delimiter=",",fmt="%1.6f")
    
