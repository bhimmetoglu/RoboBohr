# Author: Burak Himmetoglu
# Date  : 08-18-2016
# -- Project RoboBohr -- # 


import os
import shutil
import numpy as np
import commands as cm
import tarfile
from subprocess import call
from glob import glob
from Constants import *

class clusterJob:
  def __init__ (self, pathScfIn):
    # Initiate with path to scf.in files
    self.pathScfIn = pathScfIn

  def createJob(self, scheduler = "Torque", nodesize = 12, nodes=1, ppn=12, walltime="02:00:00", \
                pathPP = "./pseudos", pathPW = ".", pathWork = ".", \
                startMolecule = 0, endMolecule = 999):
    """ The method to create job submission script """

    # Number of cores
    Ncores = nodes*ppn

    # Run command
    if (scheduler == "Torque"):
      runCommand = "mpirun" + " " + "-np" + " " + str(Ncores)
    elif (scheduler == "Slurm"):
      runCommand = "ibrun"
    elif (scheduler == "None"):
      runCommand = "mpirun" + " " + "-np" + " " + str(Ncores)
    else:
      print("Stopping...")
      print("Error: Unknown scheduler")
      quit()

    # Chdir to workdir
    os.chdir(pathWork)

    # Bring pw.x
    if ( not os.path.exists("pw.x") ):
      shutil.copy2(os.path.join(pathPW,"pw.x"),"pw.x")

    # Make sure that the pseudos are in the running folder
    if ( not os.path.exists("pseudos") ):
      shutil.copytree(pathPP,"./pseudos", True) 

    # Get a list of all the input files 
    inputFiles = sorted(os.listdir(self.pathScfIn))

    # Copy the desired number of input files to present directory
    for fInp in inputFiles[startMolecule:endMolecule+1]:
      shutil.copy2(os.path.join(self.pathScfIn,fInp),fInp)
      # Create the temporary folders
      moleculeId = fInp[:-7]; tempDir = "temp." + moleculeId
      if ( not os.path.exists(tempDir) ):
        os.mkdir(tempDir)

    # Id of the Job
    jobId = "rb" + "_" + str(startMolecule) + "_" + str(endMolecule)

    # Create the submission files
    fSubName = "submit." + str(startMolecule) + "_" + str(endMolecule)
    fSub = open(fSubName,"w+")
    if ( scheduler == "Torque" ):
      fSub.write("""#!/bin/bash
#PBS -l nodes=%(nodes)i:ppn=%(ppn)i
#PBS -l walltime=%(walltime)s
#PBS -j oe

export WORK_DIR=%(pathWork)s
export NCORES=%(Ncores)i

cd $WORK_DIR

""" % vars())
    elif ( scheduler == "Slurm"):
      fSub.write("""#!/bin/bash
#SBATCH --job-name="%(jobId)s"
#SBATCH --output="rb.out"
#SBATCH --partition=compute
#SBATCH --nodes=%(nodes)i
#SBATCH --ntasks-per-node=%(ppn)i
#SBATCH --export=ALL
#SBATCH -t %(walltime)s
#SBATCH -A TG-ASC090080

export WORK_DIR=%(pathWork)s

""" % vars())

    # Add mpiruns or ibruns 
    for fInp in inputFiles[startMolecule:endMolecule+1]:
      moleculeId = fInp[:-7]
      tempDir = "temp." + moleculeId
      fOut = moleculeId + ".scf.out"
      if ( scheduler == "Torque"):
        fSub.write("""%(runCommand)s ./pw.x < %(fInp)s > %(fOut)s   
rm -r ./%(tempDir)s
mv %(fOut)s ./results
""" % vars())
      elif ( scheduler == "Slurm"): # For Slurm no need to remove temp folders, they will be purged later
        fSub.write("""%(runCommand)s ./pw.x < %(fInp)s > %(fOut)s
rm -r ./%(tempDir)s
mv %(fOut)s ./results 
""" % vars() )

    fSub.close()

    # Save in the object
    self.fSubName = fSubName

  def submitJob(self,scheduler = "Torque", pathWork = "."):
    """ The method to submit job """

    # Chdir to workdir
    os.chdir(pathWork)

    # Create a results folder if non-existent
    if ( not os.path.exists("results") ):
      os.mkdir("results")

    # Submit the job to the queue
    if (scheduler == "Torque"):
      call("qsub " + self.fSubName, shell = True)
    elif (scheduler == "Slurm"):
      call("sbatch " + self.fSubName, shell = True)


