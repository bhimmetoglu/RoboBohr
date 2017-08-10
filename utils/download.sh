#!/bin/bash
# Author: Burak Himmetoglu
# Date  : 04-13-2017
# -- AlcNet -- # 

# Download data from PubChem repository
wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound_3D/01_conf_per_cmpd/SDF/00000001_00025000.sdf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound_3D/01_conf_per_cmpd/SDF/00025001_00050000.sdf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound_3D/01_conf_per_cmpd/SDF/00050001_00075000.sdf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound_3D/01_conf_per_cmpd/SDF/00075001_00100000.sdf.gz

mkdir -o SDF
mv *sdf.gz ./SDF
