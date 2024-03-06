# Run with the command:
#   python setup.py build_ext --inplace

from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

# For HTML output of where python/C is being used.
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True


extensions = [
    Extension("*", ["*.pyx"],
        include_dirs=[numpy.get_include()]),
]
setup(
    name="flan",
    ext_modules=cythonize(extensions),
)

#setup(
#	ext_modules = cythonize(
#		["impurity_transport.pyx", "imp_cpp.pyx"],
#		annotate=True)
#)
