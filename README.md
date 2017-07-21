# RoboBohr
RoboBohr is a machine learning framework for predicting electronic structure of molecules. 
RoboBohr currently has 4 modes of operation:

1. *query*: Reads input sdf files and creates list of objects that contain types of atoms and coordinates for each entry in the sdf input. The list of these objects are then used to create input files for the pwscf code of the [Quantum Espresso](http://www.quantum-espresso.org/) package.
2. *createFeatures*: Fom the list of objects generated in the query step, generates a data matrix and saves on file.
3. *cluster*: Creates job submission files for submitting the pwscf input files in an HPC environment. Torque and Slurm scheduling systems are supported.
4. *outcomes*: Analyzes and output files generated from pwscf runs and stores relevant outcome quantities (e.g. ground state energies) and creates a log file.

### Data

The sdf data files can be downloaded from the [PubChem database](https://pubchem.ncbi.nlm.nih.gov/). Specifically, their FTP server has all the data in sdf file format. In the examples provided, the **Compound_3D/01_conf_per_cmpd/** folder has been used to generate molecular data. 

There a sereval choices of features that are available. These include, pair-distances and Coulomb matrices. The Coulomb matrix features allows random copies of the same molecule with re-shuffled indices as well as the eigenspectrum for feature matrix construction. 

### Analysis

The main part of RoboBohr is written in Python with extensions written in C (wrapped using [Cython](http://cython.org/)). The data analsysis, visualization and training of learning algorithms are performed by R scripts included in the repo.

### Article
The preprint of the article explaining RoboBohr can be accessed from [here](http://arxiv.org/abs/1609.07124).
An illustration of results from RoboBohr and alaysis performed in R can be found in the [RPubs article](https://rpubs.com/burakh/robobohr)
