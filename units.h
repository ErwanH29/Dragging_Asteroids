#ifndef UNITS_H
#define UNITS_H

#include <cmath>

const double G_si = 6.67428e-11;
const double G = 1.0;  // N-body units

class UnitSystem{
    public:
        double energy;
        double mass;
        double length;
        double time;
        double velocity;
};

extern UnitSystem units;  // Make it global

void set_units(double mass_unit_si, double length_unit_si);

// To N-body units
double mass_si_to_nb(double mass_si);
double length_si_to_nb(double length_si);
double velocity_si_to_nb(double velocity_si);
double time_si_to_nb(double time_si);
double energy_si_to_nb(double energy_so);

// To SI units
double mass_nb_to_si(double mass_nb);
double length_nb_to_si(double length_nb);
double velocity_nb_to_si(double velocity_nb);
double time_nb_to_si(double time_nb);
double energy_nb_to_si(double energy_nb);

#endif