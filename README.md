<img src="https://github.com/bhimmetoglu/RoboBohr/blob/master/roboLogo.png", width="250" height="200" />

# RoboBohr
RoboBohr is a machine learning framework for electronic structure prediction of molecules. 
RoboBohr currently has 4 modes of operation:
1. query: Reads input sdf files and creates list of objects that contain types of atoms and coordinates for each entry in the sdf input. The list of these objects are then used to create input files for the pwscf code of the [Quantum Espresso] (http://www.quantum-espresso.org/) package.
2. createFeatures: Fom the list of objects generated in the query step, generates a data matrix and saves on file.
3. cluster: Creates job submission files for submitting the pwscf input files in an HPC environment. Torque and Slurm scheduling systems are supported.
4. outcomes: Analyzes and output files generated from pwscf runs and stores relevant outcome quantities (e.g. ground state energies) and creates a log file.


