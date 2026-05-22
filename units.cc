#include "units.h"

UnitSystem units;

void set_units(double mass_unit_si, double length_unit_si){
    units.mass = mass_unit_si;
    units.length = length_unit_si;
    units.time = std::sqrt(
        units.length * units.length * units.length / (G_si * units.mass)
    );

    units.velocity = units.length / units.time;
    units.energy = units.mass * units.velocity * units.velocity;
}


double energy_si_to_nb(double energy_si){
    return energy_si / units.energy;
}


double length_si_to_nb(double length_si){
    return length_si / units.length;
}


double mass_si_to_nb(double mass_si){
    return mass_si / units.mass;
}


double time_si_to_nb(double time_si){
    return time_si / units.time;
}


double velocity_si_to_nb(double velocity_si){
    return velocity_si / units.velocity;
}


double energy_nb_to_si(double energy_nb){
    return energy_nb * units.energy;
}


double length_nb_to_si(double length_nb){
    return length_nb * units.length;
}


double mass_nb_to_si(double mass_nb){
    return mass_nb * units.mass;
}


double time_nb_to_si(double time_nb){
    return time_nb * units.time;
}


double velocity_nb_to_si(double velocity_nb){
    return velocity_nb * units.velocity;
}
