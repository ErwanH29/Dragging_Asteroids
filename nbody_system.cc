#include <cmath>
#include "nbody_system.h"
#include "units.h"

#define _TINY_ pow(2.0, -52.0)

NBodySystem::NBodySystem(){}

int NBodySystem::size(){
    return radius.size();
}

// Setters and getters for static attributes
void NBodySystem::set_eps2(double _eps2){
    eps2 = _eps2;
}

double NBodySystem::get_eps2(){
    return eps2;
}

void NBodySystem::set_mass_threshold(double _mass_threshold){
    mass_threshold = _mass_threshold;
}

double NBodySystem::get_mass_threshold(){
    return mass_threshold;
}

void NBodySystem::set_integrator_choice(int _integrator_choice){
    integrator_choice = _integrator_choice;
}

int NBodySystem::get_integrator_choice(){
    return integrator_choice;
}

// Unit conversions
void NBodySystem::set_nb_units(){
    double mass_units = 0.0;
    double length_units = 0.0;
    int nparticles = size();
    for (int i=0; i<nparticles; i++){
        mass_units += mass[i];
        for (int j=i+1; j<nparticles; j++){
            vec drij = pos[i] - pos[j];  // Using overloaded - operator
            double dr = sqrt(drij * drij);
            if (dr > length_units){
                length_units = dr;
            }
        }
    }
    set_units(mass_units, length_units);
}

void NBodySystem::convert_si_attr_to_nb(){
    int nparticles = size();
    for (int i=0; i<nparticles; i++){
        mass[i] = mass_si_to_nb(mass[i]);
        radius[i] = length_si_to_nb(radius[i]);
        pos[i] = pos[i] * (1.0 / units.length);
        vel[i] = vel[i] * (1.0 / units.velocity);
    }
}

// Add a new particle to the system
int NBodySystem::new_particle(
    double _radius,
    double _mass, 
    double x, 
    double y, 
    double z, 
    double vx, 
    double vy, 
    double vz
){
    int id = size();
    static const vec null_acc = vec(0.0, 0.0, 0.0);

    mass.push_back(_mass);
    radius.push_back(_radius);
    pos.push_back(vec(x, y, z));
    vel.push_back(vec(vx, vy, vz));
    acc.push_back(null_acc);
    pot.push_back(0.0);
    if (_mass >= mass_threshold){
        type.push_back(1);
    } else {
        type.push_back(0);
    }
    return id;
}


// Integration methods
double NBodySystem::get_energy(){
    double ke = 0.0;
    double pe = 0.0;
    for (int i=0; i<size(); i++){
        double mass_i = mass[i];
        if (mass_i <= _TINY_) continue;
        double vel2 = vel[i] * vel[i];  // Using overloaded * operator
        ke += 0.5 * mass_i * vel2;
        pe += 0.5 * mass_i * pot[i];
    }
    return ke + pe;
}

int NBodySystem::update_gravity(){
    int nparticles = size();
    static vec zero_vec = vec(0.0, 0.0, 0.0);  // Initialised at first call, and remembered throughout
    for (int i=0; i<nparticles; i++){
        acc[i] = zero_vec;
        pot[i] = 0.0;
    }

    for (int i=0; i<nparticles; i++){
        const vec& pos_i = pos[i];
        double mass_i = mass[i];
        double rad_i = radius[i];
        int type_i = type[i];
        for (int j=i+1; j<nparticles; j++){  // Current symmetric scheme can't be parallelised
            double mass_j = mass[j];
            // Test particles don't affect one another
            if ((type_i == 0) && (type[j] == 0)){
                continue;
            }

            vec rij = pos_i - pos[j];  // Using overloaded - operator
            double dr2 = rij * rij;
            double coll_rad = rad_i + radius[j];
            if (dr2 <= coll_rad * coll_rad){  // Handle collision (not yet implemented)
                continue;
            }

            dr2 += eps2;  // Softening to avoid singularities
            double dr = sqrt(dr2);
            double dr3 = dr * dr2;

            //Compute acc. and pot.
            double grav_factor = -G / dr3;
            if (mass_j>_TINY_){
                pot[i] -= G * mass_j / dr;
                double mj_factor = mass_j * grav_factor; // Order of operations can change results
                for (int k=0; k<NDIM; k++){
                    acc[i][k] += mj_factor * rij[k];
                }
            }
            if (mass_i>_TINY_){
                pot[j] -= G * mass_i / dr;
                double mi_factor = mass_i * grav_factor;
                for (int k=0; k<NDIM; k++){
                    acc[j][k] -= mi_factor * rij[k];  // drji = -drij
                }
            }
        }
    }
    return 0;
}

int NBodySystem::kick(double time_step){
    int nparticles = size();
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            vel[i][k] += time_step * acc[i][k];
        }
    }
    return 0;
}

int NBodySystem::drift(double time_step){
    int nparticles = size();
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            pos[i][k] += time_step * vel[i][k];
        }
    }
    return 0;
}

// Integration schemes
int NBodySystem::leapfrog_kdk(double time_step){
    kick(0.5 * time_step);
    drift(time_step);
    update_gravity();
    kick(0.5 * time_step);
    return 0;
}

int NBodySystem::yoshida_fourth_order(double time_step){
    // Coefficients from 1990PhLA..150..262Y
    static const double w0 = -cbrt(2.0) / (2.0 - cbrt(2.0));
    static const double w1 = 1.0 / (2.0 - cbrt(2.0));
    static const double c1 = w1 / 2.0;
    static const double c2 = (w0 + w1) / (2.0);
    static const double c3 = (w0 + w1) / (2.0);
    static const double c4 = w1 / 2.0;
    static const double d1 = w1;
    static const double d2 = w0;
    static const double d3 = w1;

    int nparticles = size();
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            pos[i][k] += c1 * time_step * vel[i][k];
        }
    }

    update_gravity();
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            vel[i][k] += d1 * time_step * acc[i][k];
        }
    }
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            pos[i][k] += c2 * time_step * vel[i][k];
        }
    }

    update_gravity();
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            vel[i][k] += d2 * time_step * acc[i][k];
        }
    }
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            pos[i][k] += c3 * time_step * vel[i][k];
        }
    }

    update_gravity();
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            vel[i][k] += d3 * time_step * acc[i][k];
        }
    }
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            pos[i][k] += c4 * time_step * vel[i][k];
        }
    }

    return 0;
}


