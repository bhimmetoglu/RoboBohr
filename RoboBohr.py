# Author: Burak Himmetoglu (burakhmmtgl@gmail.com)
# Date  : 08-18-2016
# 
# -- Project RoboBohr -- # 

import numpy as np
import sys 
import os
from Constants import *
pwd = os.path.dirname(os.path.realpath(__file__)) # Make sure to include in path 
sys.path.append(pwd); sys.path.append(bohrDir)

from classes import *
from generateQuery import *
from generateQEinputs import *
from generateData import *
from generateOut import *
from generateClusterJob import *

rBohrMode = sys.argv[1] # Read execution mode from stdin

if ("query" in rBohrMode):
  ## Initiate an instance of pubChemQuery
  query = pubChemQuery(sdfFolder = sdfPath)

  ## Read the data
  query.readData()

  ## Get listOfMolecules
  listOfMolecules = query.processQuery(elementList, valences, natMax, box, tol = 5.29)

  ## Loop over molecules
  nMolecules = len(listOfMolecules)
  f = open("AEsub.out","w+") # File that contains sum_i (n_i E_i)
  for imol in range(nMolecules):
    tempMolecule = listOfMolecules[imol]
    # Create QE inputs
    scf = QEinput(molecule = tempMolecule)
    # Write into file
    scf.scfInput(pathScfIn, foldPP = pathPP)
    # Amount to be subtracted for atomization energies:
    aEsub = scf.singleAtom()
    f.write(str(tempMolecule.sid) + " " + str(aEsub) + "\n")

  f.close()

if ("createFeatures" in rBohrMode):
  if ( featureType == "pairDistance"):
    # Store data table
    try:
      storePairFeatures(listOfMolecules)
    except Exception:
      query = pubChemQuery(sdfFolder = sdfPath)
      query.readData()
      listOfMolecules = query.processQuery(elementList, valences, natMax, box, tol = 5.29)
      storePairFeatures(listOfMolecules)

  elif ( featureType == "Coulomb"):
    ## Store data table
    try: 
      storeFeatures(listOfMolecules, eigenval = eigenval, nrandom = nrandom)
    except Exception:
      query = pubChemQuery(sdfFolder = sdfPath)
      query.readData()
      listOfMolecules = query.processQuery(elementList, valences, natMax, box, tol = 5.29)
      storeFeatures(listOfMolecules, eigenval = eigenval, nrandom = nrandom)

if ("cluster" in rBohrMode):
  ## Create input files for cluster  
  # Let us divide the whole job into chunks of 1000 molecules
  try:
    nJobs = nMolecules / 1000
  except Exception:
    query = pubChemQuery(sdfFolder = sdfPath)
    query.readData()
    listOfMolecules = query.processQuery(elementList, valences, natMax, box, tol = 5.29)
    nMolecules = len(listOfMolecules)
    nJobs = nMolecules / 1000

  # If there is roll over 
  if (nMolecules - nJobs*1000 >0):
    nJobs +=1

  # Initiate
  cluster = clusterJob(pathScfIn=pathScfIn)
  startMolecule = 0
 
  # Loop over chunks
  for job in range(nJobs):
    if (startMolecule + 999 < nMolecules):
      endMolecule = startMolecule + 999
    else:
      endMolecule = startMolecule + (nMolecules - (startMolecule+1))

    # Create inputs
    cluster.createJob(scheduler = scheduler, nodesize = nodesize, nodes=nodes, ppn=ppn, walltime=walltime, \
                     pathPP = pathPP, pathPW = pathPW, pathWork = pathWork, \
                     startMolecule = startMolecule, endMolecule = endMolecule) 

    # Reset startMolecule
    startMolecule = endMolecule+1

if ("outcomes" in rBohrMode):
  # List of outputs
  listScfOut = os.listdir(pathOut)
  nMolecules = len(listScfOut)

  ## Process scf.out queries
  query = queryScfOut(pathOut = pathOut)
  # Check outputs
  query.checkOut(startMolecule = 0, endMolecule = nMolecules)
  # Write Logs
  query.writeLogs() 

