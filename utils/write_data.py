# Author: Burak Himmetoglu
# Date in   : 04-13-2017
# Date curr : 08-10-2017
# -- AlcNet -- # 

from Constants import *
from generateQuery import *
import os

# Get a list of the sdf files
listSDF = os.listdir("./SDF/")
listSDF.sort() 

# Maximum number of atoms allowed in a given molecule
natMax = 10000

# Loop over SDF files and save
for fil in listSDF:
  print("Processing "+ fil)

  # Read the sdf file
  path = os.path.join("./SDF", fil)
  query = pubChemQuery(sdfFile=path)

  print("Parsing file " + fil + "...\n")
  query.readData()

  # Process query and save
  listOfMolecules = query.processQuery(periodic_table, return_dictionary=True, natMax=natMax)

  # Store JSON
  filename = 'pubChem_p_' + fil[:-7] + '.json'
  storeJSON(listOfMolecules, filename)
