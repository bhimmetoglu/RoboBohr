# Author: Burak Himmetoglu
# Date  : 08-18-2016
# -- Project RoboBohr -- # 

import numpy as np
import os
import re
from collections import Counter
from classes import *
from Constants import *

class QEinput:
  def __init__ (self, molecule):
    # Initiate with an instance of class molecule
    self.molecule = molecule

  def singleAtom(self):
    """ Amount to be subtracted using single atom energies """
    elements = self.molecule.names # Names of elements
    natByType = []; AE = []
    for ie in range(len(elementList)):
      listOfIndices = [ i for i,x in enumerate(elements) if x == elementList[ie]]
      natByType.append(len(listOfIndices))
      AE.append(singleAtomEnergies[ie])

    return np.dot(np.array(natByType),np.array(AE))

  def scfInput(self, pathScfIn, foldPP = ".", box = 30.0):
    # Make sure there is an outDir path
    if (not os.path.exists(pathScfIn)):
      os.mkdir(pathScfIn)

    # Get the molecular structure 
    elements = self.molecule.names # Names of elements
    pos = np.array(self.molecule.positions, dtype = float) # Atomic positions  
    typ = Counter(elements).keys() # Atomic types, i.e. unique names
    ntype = len(typ) # Number of distinct types  
    natoms = len(elements) # Total number of atoms in the molecule

    # Match masses, valences, and pseudos
    M = []; PP = []
    for it in range(ntype):
      elementIndex = elementList.index(typ[it])
      M.append(masses[elementIndex]) # Mass of the corresponding type
      PP.append(pseudos[elementIndex]) # Pseudo of the corresponding type
    
    # Find number of atoms per type (in the order of elementList), their valences and single atom energies
    natByType = []; val = []
    for ie in range(len(elementList)):
      listOfIndices = [ i for i,x in enumerate(elements) if x == elementList[ie]]
      natByType.append(len(listOfIndices))
      val.append(valences[ie])

    # Compute the number of states (nbnd)
    nbnd = np.dot(np.array(natByType),np.array(val))
    nbnd = nbnd/2 + 1

    # Standard input file 
    moleculeId = "000000000" + str(self.molecule.sid); moleculeId = moleculeId[-9:] # Consistent 9 digit naming
    filnam = os.path.join(pathScfIn,moleculeId + ".scf.in")
    tempdir = "temp." + moleculeId
    f=open(filnam,"w+")

    # Some parameters
    celldm = box
    ecut = 30.0; ecutrho = 10*ecut

    # Write CONTROL, SYSTEM and ELECTRONS cards
    f.write("""&CONTROL
 calculation = 'scf'
 outdir = '%(tempdir)s'
 prefix = '%(moleculeId)s'
 pseudo_dir = '%(foldPP)s'
/
&SYSTEM
 ibrav = 1,
 nosym = .true.,
 celldm(1) = %(box)2.1f,
 nat = %(natoms)i,
 ntyp = %(ntype)i,
 ecutwfc = %(ecut)2.1f,
 ecutrho = %(ecutrho)2.1f,
 nbnd    = %(nbnd)i,
/
&ELECTRONS
 conv_thr = 1.0d-6,
 mixing_beta = 0.3,
/
""" % vars())

    # Write ATOMIC_SPECIES card
    f.write("ATOMIC_SPECIES\n")
    for it in range(ntype):
      text = typ[it] + " " + str(M[it]) + " " + PP[it] + "\n"
      f.write(text)

    # Write ATOMIC_POSITIONS card
    f.write("ATOMIC_POSITIONS angstrom \n")
    for iat in range(natoms):
      text = elements[iat] + " " + str("%.8f" % pos[iat,0]) + " " + str("%.8f" % pos[iat,1]) + " " + str("%.8f" % pos[iat,2]) + "\n"
      f.write(text)

    # Write K_POINTS card
    f.write("K_POINTS gamma")

    # Close
    f.close()

