# distutils: language=c++

# Copying the example found here:
# https://cython.readthedocs.io/en/stable/src/userguide/wrapping_CPlusPlus.html#wrapping-cplusplus

from libcpp.vector cimport vector

cdef extern from "Impurity.cpp":
    pass

# Declare the class with cdef
cdef extern from "Impurity.h" namespace "impurities":
    cdef cppclass Impurity:
        Impurity() except +
        Impurity(int, float, float, float, float, int, int, float) except +
        int fstart, imp_atom_num, charge;
        float mass, x, y, z, vx, vy, vz, tz, t, xstart, zstart, weight;
        vector[float] xhist, yhist, zhist, vxhist, vyhist, vzhist;
