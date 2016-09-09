from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext

# Build: python setup.py build_ext -i
setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("pairFeatures",
                 sources=["_pairFeatures.pyx", "pairFeatures.c"],
                 extra_compile_args=['-std=gnu99'],
                 include_dirs=[numpy.get_include()])],
)

