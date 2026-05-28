#ifndef UNITS_H
#define UNITS_H

#include <cmath>

inline double _TINY_ = pow(2.0, -52.0);
constexpr double pi2 = M_PI * M_PI;
constexpr double G_si = 6.67428e-11;
constexpr double G = 1.0;  // N-body units

class UnitSystem{
    public:
        double energy;
        double mass;
        double length;
        double time;
        double velocity;

        void set_units(double mass_unit_si, double length_unit_si);

        double mass_si_to_nb(double mass_si);
        double length_si_to_nb(double length_si);
        double velocity_si_to_nb(double velocity_si);
        double time_si_to_nb(double time_si);
        double energy_si_to_nb(double energy_si);

        double mass_nb_to_si(double mass_nb);
        double length_nb_to_si(double length_nb);
        double velocity_nb_to_si(double velocity_nb);
        double time_nb_to_si(double time_nb);
        double energy_nb_to_si(double energy_nb);
};

extern UnitSystem units;  // Make it global, one instance exists in .cc


#endif