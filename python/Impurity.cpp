#include <iostream>
#include "Impurity.h"

// Copying the example found here:
// https://cython.readthedocs.io/en/stable/src/userguide/wrapping_CPlusPlus.html#wrapping-cplusplus
namespace impurities {

    // Default constructor
    Impurity::Impurity () {}

    // Overloaded constructor
    Impurity::Impurity (int imp_atom_num, float mass, float x, float y, 
                float z, int charge, int fstart, float weight) {
        this->imp_atom_num = imp_atom_num;
        this->charge = charge;
        this->mass = mass;
        this->x=x;
        this->y=y;
        this->z=z;
        this->fstart=fstart;
        this->zstart=z;
        this->xstart=x;
        this->weight=weight;
        this->vx=0;
        this->vy=0;
        this->vz=0;
        this->tz=0;
        this->t=0;
    
    }

    // Destructor
    Impurity::~Impurity () {}
}
