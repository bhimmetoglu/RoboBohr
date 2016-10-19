# Author: Burak Himmetoglu
# Date  : 08-18-2016
# -- Project RoboBohr -- # 

import numpy as np

## Example sdf 
#ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound_3D/01_conf_per_cmpd/SDF/00050001_00075000.sdf.gz

## Some constants
auToA = 0.529 # au is 0.529 A

## Paths to specific directories
bohrDir = "/home/burak/NielsBohr/RoboBohr/src" # Path to the location of RoboBohr
sdfPath = "/home/burak/NielsBohr/sdf" # Path to the folder where sdf files are located

# Data input/output
pathScfIn = "/home/burak/NielsBohr/scf" # Location where scf.in files will be stored
pathOut = "/home/burak/NielsBohr/results" # Location where scf.out files will be stores
pathData = "/home/burak/NielsBohr/data" # Location where data will be written

## Parameters for RoboBohr 
# Here is an example for studying CHNOPS
elementList = [ "H", "C", "N", "O", "P", "S" ]
valences = [1, 4, 5, 6, 5, 6]
masses = [1.0, 12.0, 14.0, 15.999, 30.97, 32.06]
ZList = [1, 6, 7, 8, 15, 16]
pseudos = ["H.pbe-rrkjus_psl.0.1.UPF", \
           "C.pbe-n-rrkjus_psl.0.1.UPF", \
           "N.pbe-n-rrkjus_psl.0.1.UPF", \
           "O.pbe-n-rrkjus_psl.0.1.UPF", \
           "P.pbe-n-rrkjus_psl.0.1.UPF", \
           "S.pbe-n-rrkjus_psl.0.1.UPF" ]
singleAtomEnergies = [  -0.917798,\
                        -11.275352,\
                        -19.577858,\
                        -33.171237,\
                        -15.12336,\
                        -22.522057] # Obtained from ld1, from PP generation output

## Cluster settings
nMolecules = 16273 
scheduler = "Torque"; nodesize = 12; nodes=1; ppn=12; walltime="02:00:00"

# Location of certain files
pathPP = "/home/burak/pseudos" # Location of pseudopotentials 
pathPW = "/home/burak/espresso/bin" # Location of pw.x  
pathWork = "/home/burak/tempDir" # Location where calculations will be performed (could be a temp location in the cluster) 


## What type of features do you want?
featureType = "Coulomb" # The other option is pairDistance

## Some parameters for generating data matrices
natMax = 50; box = 30;
eigenval = True; nrandom = 0 # Useful for Coulomb
