# Author: Burak Himmetoglu
# Date  : 08-18-2016
# -- Project RoboBohr -- # 

import numpy as np
import re
import os
import shutil
import tarfile
from subprocess import call
from glob import glob
from datetime import datetime
from Constants import *

class queryScfOut:
  """ A class for processing scf.out files """
  def __init__ (self,pathOut):
    # Initiate with pathOut, which contains the scf.out files
    self.pathOut = pathOut

  def checkOut(self, startMolecule, endMolecule):
    """ startMolecule : Id of first molecule in pathOut
        endMolecule   : Id of last molecule in pathOut
    """
    # Initiate
    self.startMolecule = startMolecule; self.endMolecule = endMolecule

    # Chdir to pathOut
    os.chdir(self.pathOut)

    # List of input/output files
    listFilesIn = sorted(glob("*.scf.in"))  
    listFilesOut = sorted(glob("*.scf.out"))  

    def stripId(x):
      # Removes .scf.in and .scf.out. x has to be a list or array
      return x[:9]

    # Check whether the number of files are equal. 
    # If there are unfinished runs, store them in a log
    if ( len(listFilesIn) != len(listFilesOut) ):
      # If there are scf.in's that are not computed
      idIn = map(stripId, listFilesIn)
      idOut = map(stripId, listFilesOut)
      diff = list(set(idIn) - set(idOut))
      f = open("unfinished.log","w+")
      for ids in diff:
        fil = ids + ".scf.in"
        f.write(fil + "\n")
      f.close()
 
    # Save in object    
    self.listFilesIn = listFilesIn; self.listFilesOut = listFilesOut

  def writeLogs(self):
    # Chdir to pathOut
    os.chdir(self.pathOut)
    # Create a log file to keep track of which computations are performed
    logFile = "run.log" + "." + str(self.startMolecule) + "_" + str(self.endMolecule)
    # Create an output data file to store energies
    dataFile = "out.dat" + "." + str(self.startMolecule) + "_" + str(self.endMolecule)
    # Open files
    fLog = open(logFile,"w+"); fDat = open(dataFile,"w+")
    # Write the date/time and start/end ids
    fLog.write(str(datetime.today()) + "\n") # Write the date and time on top
    fLog.write("start: " + " " +  str(self.startMolecule) + " " + " " + "end: " + " " + str(self.endMolecule) + "\n")
    fDat.write(str(datetime.today()) + "\n") # Write the date and time on top
    fDat.write("start: " + " " +  str(self.startMolecule) + " " + " " + "end: " + " " + str(self.endMolecule) + "\n")

    # Loop over output files
    for filOut in self.listFilesOut:
      moleculeId = filOut[:9]
      #filOut = moleculeId + ".scf.out"
      # Read the output
      f = open(filOut,"r")
      tempOut = f.readlines()
      f.close()
      # Get output energy
      lConv = False 
      for line in tempOut:
        line = line.rstrip()
        if re.search("!", line):
          fLog.write("""%(filOut)s:  done \n""" % vars())
          lConv = True
          E = line.split()[4] # Ground state energy
          break
      # No convergence:
      if (lConv == False):
        fLog.write("""%(filOut)s:  no_convergence \n""" % vars())
        
      # Write into fDat
      if (lConv == True):
        fDat.write(moleculeId + " " + E + "\n")
      else:
        fDat.write(moleculeId + " " + "NA" + "\n")
      
    fLog.close()
    fDat.close()

  def archiveFiles(self, destination):
    # Archive files 
    filnam = "results" + "." + str(self.startMolecule) + "_" + str(self.endMolecule) + ".tar.gz"
    tar = tarfile.open(filnam, "w:gz")
    for fil in os.listdir(self.pathOut):
      tar.add(fil)
    tar.close()
    shutil.copy2(filnam,destination) # Copy to desired destination
