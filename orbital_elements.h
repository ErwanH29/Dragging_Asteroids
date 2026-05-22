#ifndef ORBITAL_ELEMENTS
#define ORBITAL_ELEMENTS

#include <cmath>

#include "units.h"
#include "vec_3d.h"

vec cross_product(vec a, vec b);

const double get_orbital_energy(
    double mass1, double mass2, 
    vec pos1, vec pos2, 
    vec vec1, vec vec2
);

const double get_sma(
    double mass1, double mass2, 
    vec pos1, vec pos2, 
    vec vec1, vec vec2
);

const double get_ecc(
    double mass1, double mass2, 
    vec pos1, vec pos2, 
    vec vec1, vec vec2
);


#endif