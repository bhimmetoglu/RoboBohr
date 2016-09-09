""" Wrapper for pairFeatures.c """

# cimport the Cython declarations for numpy
cimport numpy as np

# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example, but good practice)
np.import_array()

# cdefine the signature of our c functions
cdef extern from "pairFeatures.h":
  void pairFeatures(int nSpecies, int natoms, int * elements, double * pos, double * featMat)

# Wrapper Code, with numpy type annotations
def fun_pairFeatures(int nSpecies, int natoms,
                     np.ndarray[int, ndim=1, mode="c"] elements not None,
                     np.ndarray[double, ndim=2, mode="c"] pos not None,
                     np.ndarray[double, ndim=2, mode="c"] featMat not None):
    pairFeatures(nSpecies, natoms,
                 <int*> np.PyArray_DATA(elements),
                 <double*> np.PyArray_DATA(pos),
                 <double*> np.PyArray_DATA(featMat))

