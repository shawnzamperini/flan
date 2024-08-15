#ifndef IMPURITY_H
#define IMPURITY_H

#include <vector>


// Copying the example found here:
// https://cython.readthedocs.io/en/stable/src/userguide/wrapping_CPlusPlus.html#wrapping-cplusplus
namespace impurities {
    class Impurity {
        public:
	    int fstart, imp_atom_num, charge;
            float mass, x, y, z, vx, vy, vz, tz, t, xstart, zstart, weight;
            std::vector<float> xhist, yhist, zhist, vxhist, vyhist, vzhist;
            Impurity();
            Impurity(int imp_atom_num, float mass, float x, float y, 
				float z, int charge, int fstart, float weight);
            ~Impurity();
    };
}

#endif
