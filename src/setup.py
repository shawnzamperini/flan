# Run with the command:
#   python setup.py build_ext --inplace

from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy
import sys
import os
from distutils import sysconfig

# For HTML output of where python/C is being used.
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

# Want to use the actual GNU compiler. 
os.environ["CC"] = "g++-13"
os.environ["CXX"] = "g++-13"


if sys.platform.startswith("win"):
    openmp_arg = '/openmp'
else:
    openmp_arg = '-fopenmp'



extensions = [
    Extension("*", ["*.pyx"],
        include_dirs=[numpy.get_include()]
        #extra_compile_args=[openmp_arg],
        #extra_link_args=[openmp_arg]
    ),
]



setup(
    name="flan",
    ext_modules=cythonize(extensions),
)

# me make damn sure, that disutils does not mess with our
# build process


