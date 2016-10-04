# Author: Burak Himmetoglu
# Date  : 09-29-2016
# -- Project RoboBohr -- # 

import numpy as np
import os
import re
import random

## Some parameters
nDim = 23 # Max dimension of the Coulomb matrix

# Read the coulomb.txt file
f = open("coulomb.txt","r")
fParse = f.readlines()
f.close()

# Read properties.txt
f = open("properties.txt","r")
pParse = f.readlines()
f.close()

## Outcomes

ind = 0; out = []
## Read Properties
for line in pParse:
  temp = line.split()
  out.append(temp) 
  ind += 1

nMolecules = ind 
outcomes = np.array(out, dtype = float) # Save as array

## Coulomb

# Loop over all the lines
ind = 0; indMol = 0
temp = []
CM = np.zeros((nMolecules,nDim*nDim)) # Initiate Coulomb Matrix

for line in fParse:
  # If break is reached pass to next molecule
  if ( (nDim+1)*indMol + nDim == ind):
    CM[indMol,:] = map(float,temp) # Save the matrix (unrolled)
    temp = [] # Reset temp
    ind += 1
    indMol += 1
    continue 

  # Append matrix
  temp += line.split()

  # Next line
  ind += 1

if (indMol != nMolecules):
  print "Wrong number of molecules!"
  quit()

# Now that the Coulomb matrix for each molecule is read, we look for the upper triangular part
nTot = nDim * (nDim+1)/2 # Size of upper triangular form
iu = np.triu_indices(nDim) # Indices of upper triangular part

## Final data

# Initialize
X = np.zeros((nMolecules,nTot)) # Final data matrix (Coulomb matrix with ordered)
XL = np.zeros((nMolecules,nDim)) # Final data matrix (eigenvalues)

random.seed(101) 
for indMol in range(nMolecules):
  cmTemp = np.reshape(CM[indMol,:],(nDim,nDim)) # Refold into matrix
  # Check that it is a square matrix
  if (np.any(cmTemp-np.transpose(cmTemp)) > 1e-6):
    print "Wrong symmetry in Coulomb matrix for molecule:", indMol
    print cmTemp
    quit()

  cn = np.apply_along_axis(np.linalg.norm,1,cmTemp) # Column norms 

  ## Simply ordered matrices (no random reordering)
  Ctemp = -cn; P = Ctemp.argsort() # Get the permutation
  Xtemp = (cmTemp[:,P])[P,:] # Reordered matrix
  X[indMol,:] = Xtemp[iu] # Upper triangular part

  ## Eigenvalues
  w,v = np.linalg.eig(cmTemp)
  XL[indMol,:] = np.sort(w.real)

## Save csv files
filNam = "coulombX" + ".csv"
np.savetxt(filNam,X,delimiter=",",fmt='%1.6f')

filNam = "coulombL" + ".csv"
np.savetxt(filNam,XL,delimiter=",",fmt='%1.6f')

filNam = "outcomes" + ".csv"
np.savetxt(filNam,outcomes,delimiter=",",fmt='%1.6f') 
