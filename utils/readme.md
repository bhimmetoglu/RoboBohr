# Predict molecular properties

## Usage
Download data from PubChem ftp site using `download.sh`. Change the number of files you would like to download.

Use `write_data.py` to store in JSON format. Then, use the json library to read the data, for example:

```python
with open("pubChem_p_00025001_00050000.json") as f:
    data = json.load(f)
```

## Context
This dataset contains molecular properties scraped from the [PubChem database](http://pubchem.ncbi.nlm.nih.gov).
Each file contains properties for 25,000 molecules, made up of the elements H, C, N, O, F, Si, P, S, Cl, Br, and I.

## Data Description

Each JSON file contains a list of molecular data. An example molecule is given below:

```javascript
{'En': 37.801,
 'atoms': [{'type': 'O', 'xyz': [0.3387, 0.9262, 0.46]},
  {'type': 'O', 'xyz': [3.4786, -1.7069, -0.3119]},
  {'type': 'O', 'xyz': [1.8428, -1.4073, 1.2523]},
  {'type': 'O', 'xyz': [0.4166, 2.5213, -1.2091]},
  {'type': 'N', 'xyz': [-2.2359, -0.7251, 0.027]},
  {'type': 'C', 'xyz': [-0.7783, -1.1579, 0.0914]},
  {'type': 'C', 'xyz': [0.1368, -0.0961, -0.5161]},
  {'type': 'C', 'xyz': [-3.1119, -1.7972, 0.659]},
  {'type': 'C', 'xyz': [-2.4103, 0.5837, 0.784]},
  {'type': 'C', 'xyz': [-2.6433, -0.5289, -1.426]},
  {'type': 'C', 'xyz': [1.4879, -0.6438, -0.9795]},
  {'type': 'C', 'xyz': [2.3478, -1.3163, 0.1002]},
  {'type': 'C', 'xyz': [0.4627, 2.1935, -0.0312]},
  {'type': 'C', 'xyz': [0.6678, 3.1549, 1.1001]},
  {'type': 'H', 'xyz': [-0.7073, -2.1051, -0.4563]},
  {'type': 'H', 'xyz': [-0.5669, -1.3392, 1.1503]},
  {'type': 'H', 'xyz': [-0.3089, 0.3239, -1.4193]},
  {'type': 'H', 'xyz': [-2.9705, -2.7295, 0.1044]},
  {'type': 'H', 'xyz': [-2.8083, -1.921, 1.7028]},
  {'type': 'H', 'xyz': [-4.1563, -1.4762, 0.6031]},
  {'type': 'H', 'xyz': [-2.0398, 1.417, 0.1863]},
  {'type': 'H', 'xyz': [-3.4837, 0.7378, 0.9384]},
  {'type': 'H', 'xyz': [-1.9129, 0.5071, 1.7551]},
  {'type': 'H', 'xyz': [-2.245, 0.4089, -1.819]},
  {'type': 'H', 'xyz': [-2.3, -1.3879, -2.01]},
  {'type': 'H', 'xyz': [-3.7365, -0.4723, -1.463]},
  {'type': 'H', 'xyz': [1.3299, -1.3744, -1.7823]},
  {'type': 'H', 'xyz': [2.09, 0.1756, -1.3923]},
  {'type': 'H', 'xyz': [-0.1953, 3.128, 1.7699]},
  {'type': 'H', 'xyz': [0.7681, 4.1684, 0.7012]},
  {'type': 'H', 'xyz': [1.5832, 2.901, 1.6404]}],
 'id': 1,
 'shapeM': [259.66,
  4.28,
  3.04,
  1.21,
  1.75,
  2.55,
  0.16,
  -3.13,
  -0.22,
  -2.18,
  -0.56,
  0.21,
  0.17,
  0.09]}

```

* En: This field is the molecular energy calculated using a force-field method. See references [1,2] for details. This is the target variable which is being predicted.
* atoms: This field contains the name of the element and the position (x,y,z coordinates) and needs to be used for feature engineering.
* id : This field is the PubChem Id 
* shapeM : This field contains the shape multipoles and can be used for feature engineering. For definition of shape multipoles, see reference [3].

Notice that each molecule contains different number and types of atoms, so it is challenging to come up with features that can describe every molecule in 
a unique way. There are several approaches taken in the literature (see the references), one of which is to use the Coulomb Matrix for a given molecule
defined by

$$
C_{IJ} = \frac{Z_I\, Z_J}{\vert R_I - R_J \vert}, \quad  ({\rm I \neq J}) \qquad
C_{IJ} = Z_I^{2.4}, \quad (I=J)
$$

where $Z_I$ are atomic numbers (can be looked up from the periodic table for each element), and ${\vert R_I - R_J \vert}$ is the distance between two atoms I and J. 

## Inspiration
Simulations of molecular properties are computationally expensive. The purpose of this project is to use machine learning methods to come up with a model that can predict molecular properties from a database. In the PubChem database, there are around 100,000,000 molecules. It could take years to do simulations on all 
of these molecules, however machine learning can be used to predict their properties much faster. As a result, this could open up many possibilities in computational design and discovery of molecules, compounds and new drugs.

This is a regression problem so mean squared error is minimized during training.

## References
[1] Halgren TA. Merck Molecular Force Field: I. Basis, Form, Scope,
    Parameterization and Performance of MMFF94.  J. Comp. Chem.
    1996;17:490-519.
    
[2] Halgren TA. Merck Molecular Force Field: VI. MMFF94s Option for
    Energy Minimization Studies.  J. Comp. Chem. 1999;20:720-729.
    
[3] Kim, Sunghwan, Evan E Bolton, and Stephen H Bryant. “PubChem3D: Shape Compatibility Filtering Using Molecular Shape Quadrupoles.” J. Cheminf. 3 (2011): 25.

[4] Himmetoglu B.: Tree based machine learning framework for predicting ground state energies of molecules, J. Chem. Phys 145, 134101 (2016)

[5] Rupp M., Ramakrishnan R., von Lilienfeld OA.: Machine Learning for Quantum Mechanical Properties of Atoms in Molecules, J. Phys. Chem. Lett. , 6(16): 3309–3313 (2015)
	
[6] Montavon G., Rupp M., Gobre V., Vazquez-Mayagoitia A., Hansen K., Tkatchenko A., Müller K-R., von Lilienfeld OA.: Machine learning of molecular electronic properties in chemical compound space, New J. Phys., 15(9): 095003 (2013)

[7] Hansen K., Montavon G., Biegler F., Fazli S., Rupp M., Scheffler M., von Lilienfeld OA., Tkatchenko A., Müller K-P.: Assessment and Validation of Machine Learning Methods for Predicting Molecular Atomization Energies, J. Chem. Theory Comput. , 9(8): 3543–3556 (2013)
