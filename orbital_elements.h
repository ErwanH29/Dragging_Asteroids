#ifndef ORBITAL_ELEMENTS
#define ORBITAL_ELEMENTS

#include <cmath>

#include "units.h"
#include "vec_3d.h"


double get_orbital_energy(
    double mass1, 
    double mass2, 
    const vec& pos1, const vec& pos2, 
    const vec& vec1, const vec& vec2
);

double get_sma(
    double mass1,
    double mass2, 
    const vec& pos1, const vec& pos2, 
    const vec& vec1, const vec& vec2
);

double get_ecc(
    double mass1,
    double mass2, 
    const vec& pos1, const vec& pos2, 
    const vec& vec1, const vec& vec2
);

double get_orbital_period(
    double mass1,
    double mass2,
    double sma
);


#endif