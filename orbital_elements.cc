#include <cmath>
#include <math.h>

#include "units.h"
#include "vec_3d.h"

// Only seen locally, prevent overwritting if another function with same name is made
static vec cross_product(const vec& a, const vec& b){
    return vec(
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    );
}

double get_orbital_period(double mass1, double mass2, double sma){
    const double sma3 = sma * sma * sma;
    const double pi2 = M_PI * M_PI;
    double period_sq = 4.0 * sma3 * pi2  / (G * (mass1 + mass2));
    return sqrt(period_sq);
}

double get_orbital_energy(
    double mass1, double mass2, 
    const vec& pos1, const vec& pos2, 
    const vec& vec1, const vec& vec2
){
    vec vij = vec1 - vec2;
    vec rij = pos1 - pos2;
    double v2 = vij * vij;
    double r = sqrt(rij * rij);
    return 0.5 * v2 - G * (mass1 + mass2) / r;
}

double get_sma(
    double mass1, double mass2, 
    const vec& pos1, const vec& pos2, 
    const vec& vec1, const vec& vec2
){
    double mu = G * (mass1 + mass2);
    double energy = get_orbital_energy(
        mass1, mass2, pos1, pos2, vec1, vec2
    );
    return -mu / (2.0 * energy);
}

double get_ecc(
    double mass1, double mass2, 
    const vec& pos1, const vec& pos2, 
    const vec& vec1, const vec& vec2
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
}
