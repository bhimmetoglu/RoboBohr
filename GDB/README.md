## Analysis of the GDB data

In this folder, we provide the scripts for analyzing the GDB dataset of molecules used in Montavon et al. New J. Phys. 15 095003 (2013).
The data can be downloaded from the New Journal of Physics website as the [supplemental data](http://iopscience.iop.org/1367-2630/15/9/095003/media).

The script **readCoulomb.py** reads the Coulomb matrices and the outcomes of DFT calculations and stores them in CSV format. 
The R scripts perform the training of the boosted regression tree algorithm. **analyzeNJP_X.R** use the ordered Coulomb matrices, while **analyzeNJP_L.R** use the eigenspectrum. From the scripts, we observe two results:

1. The boosted regression tree algorithm perform as good as convolution neural networks for this dataset (compare RMSE results with the paper. Make sure to convert eV units in the paper to kcal/mol)
2. Compared to the PubChem data we used, the GDB dataset has much lower variability, leading to smaller RMSE values.
