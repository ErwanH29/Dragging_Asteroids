#ifndef ORBITAL_ELEMENTS
#define ORBITAL_ELEMENTS

#include <cmath>

#include "units.h"
#include "vec_3d.h"

vec cross_product(vec a, vec b){
    return vec(
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    );
};

const double get_orbital_energy(
    double mass1, double mass2, 
    vec pos1, vec pos2, 
    vec vec1, vec vec2
){
    vec vij = vec1 - vec2;  // Relative velocity
    vec rij = pos1 - pos2;
    double v2 = vij * vij;
    double r = sqrt(rij * rij);
    return 0.5 * v2 - G * (mass1 + mass2) / r;
};

const double get_sma(
    double mass1, double mass2, 
    vec pos1, vec pos2, 
    vec vec1, vec vec2
){
    double mu = -G * (mass1 + mass2);
    double energy = get_orbital_energy(
        mass1, mass2, pos1, pos2, vec1, vec2
    );
    return -mu / (2 * energy);
};

const double get_ecc(
    double mass1, double mass2, 
    vec pos1, vec pos2, 
    vec vec1, vec vec2
){
    vec rij = pos1 - pos2;
    vec vij = vec1 - vec2;
    vec h = cross_product(rij, vij);  // Specific angular momentum
    double h2 = h * h;
    double eps = get_orbital_energy(
        mass1, mass2, pos1, pos2, vec1, vec2
    );
    double mu = G * (mass1 + mass2);
    return sqrt(1.0 + 2.0 * eps * h2 / (mu * mu));
};


#endif