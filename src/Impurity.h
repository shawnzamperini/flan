#ifndef IMPURITY_H
#define IMPURITY_H

#include <vector>


// Copying the example found here:
// https://cython.readthedocs.io/en/stable/src/userguide/wrapping_CPlusPlus.html#wrapping-cplusplus
namespace impurities {
    class Impurity {
        public:
			int fbirth, imp_atom_num;
            float charge, mass, x, y, z, vx, vy, vz, t;
            std::vector<float> xhist, yhist, zhist, vxhist, vyhist, vzhist;
            Impurity();
            Impurity(int imp_atom_num, float mass, float x, float y, 
				float z, float charge, int fbirth);
            ~Impurity();
    };
}

#endif
