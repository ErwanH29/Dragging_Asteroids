/* 
Symplectic integrator with 1st, 2nd, 3rd and 4th
order schemes. Symplectic integrators preserve
energy error over secular times, however across
an orbital period some error does arise. 

This error can be reduced by making the internal 
time-step an integer fraction of the orbital period. 
This error can be reduced further by making it an 
integer fraction of the numerical orbital period, 
and further reduced by making this a power of two 
integer fraction of orbital period (numerical round off).
*/

#include <chrono>
#include <cmath>
#include <vector>
#include "units.h"
#include "vec_3d.h"

#define _TINY_ pow(2.0, -52.0)

template <typename T>
using dyn_arr = std::vector<T>;

const int NDIM = 3;
static const double eps2 = 0.0;
static dyn_arr<double> mass, pot, radius;
static dyn_arr<vec> acc, pos, vel;


void set_nb_units(){
    double mass_units = 0.0;
    double length_units = 0.0;
    int nparticles = radius.size();
    for (int i=0; i<nparticles; i++){
        mass_units += mass[i];
        for (int j=i+1; j<nparticles; j++){
            vec pos_i = pos[i];
            vec pos_j = pos[j];
            vec drij = pos_i - pos_j;  // Using overloaded - operator
            double dr = sqrt(drij * drij);
            if (dr > length_units){
                length_units = dr;
            }
        }
    }
    set_units(mass_units, length_units);
};


void convert_particles_to_nb(){
    int nparticles = radius.size();
    for (int i=0; i<nparticles; i++){
        mass[i] = mass_si_to_nb(mass[i]);
        radius[i] = length_si_to_nb(radius[i]);
        pos[i] = pos[i] * (1.0 / units.length);
        vel[i] = vel[i] * (1.0 / units.velocity);
    }
}


int new_particle( 
    double _radius,
    double _mass, 
    double x, double y, double z, 
    double vx, double vy, double vz
){
    int id = radius.size();
    mass.push_back(_mass);
    radius.push_back(_radius);
    pos.push_back(vec(x, y, z));
    vel.push_back(vec(vx, vy, vz));
    acc.push_back(vec(0.0, 0.0, 0.0));
    pot.push_back(0.0);
    return id;
};


double get_energy(){
    int Nparticles = radius.size();
    double ke = 0.0;
    double pe = 0.0;
    for (int i=0; i<Nparticles; i++){
        double mass_i = mass[i];
        if (mass_i <= _TINY_) continue;
        double vel2 = vel[i] * vel[i];  // Using overloaded * operator
        ke += 0.5 * mass_i * vel2;
        pe -= 0.5 * G * mass_i * pot[i];
    }
    return ke + pe;
};


int get_acc_pot_coll(){
    int Nparticles = radius.size();
    for (int i=0; i<Nparticles; i++){
        acc[i] = vec(0.0, 0.0, 0.0);
        pot[i] = 0.0;
    }

    for (int i=0; i<Nparticles; i++){
        vec pos_i = pos[i];
        vec vel_i = vel[i];
        double mass_i = mass[i];
        double rad_i = radius[i];
        for (int j=i+1; j<Nparticles; j++){
            double mass_j = mass[j];
            // Test particles don't affect one another
            if ((mass_j <= _TINY_) && (mass_i <= _TINY_)){
                continue;
            }

            vec pos_j = pos[j];
            vec vel_j = vel[j];
            vec rij = pos_i - pos_j;  // Using overloaded - operator
            vec vij = vel_i - vel_j;
            double dr2 = rij * rij;
            double coll_rad = rad_i + radius[j];
            if (dr2 <= coll_rad * coll_rad){
                // Handle collision (not yet implemented)
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
};


int leapfrog_kdk(double time_step){
    int Nparticles = radius.size();

    // 1/2 kick step
    for (int i=0; i<Nparticles; i++){
        for (int k=0; k<NDIM; k++){
            vel[i][k] += 0.5 * time_step * acc[i][k];
        }
    }

    // drift step
    for (int i=0; i<Nparticles; i++){
        for (int k=0; k<NDIM; k++){
            pos[i][k] += time_step * vel[i][k];
        }
    }

    // 2/2 kick step
    get_acc_pot_coll();  // Update acc.
    for (int i=0; i<Nparticles; i++){
        for (int k=0; k<NDIM; k++){
            vel[i][k] += 0.5 * time_step * acc[i][k];
        }
    }

    return 0;
}


int main(){
    double M_si = 2e30;                // Mass of the sun
    double r_si = 1.5e11;              // 1 AU
    double vorb_si = sqrt(G_si * M_si / r_si);  // Orbital velocity

    int p1 = new_particle(
        1.0, M_si, 
        0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0
    );
    int p2 = new_particle(
        1.0, M_si, 
        -r_si, 0.0, 0.0, 
        0.0, -0.99 * vorb_si, 0.0
    );

    set_nb_units();
    convert_particles_to_nb();

    double time = 0.0;
    const double end_time_si = 10000.0 * 3600 * 24 * 365.25;  // 100 years
    const double end_time_nb = time_si_to_nb(end_time_si);
    const double time_step = end_time_nb / pow(2.0, 20.0);  // 32k steps per orbit

    get_acc_pot_coll();
    const double initial_energy = get_energy();

    std::cout << "MASS_UNIT: " << units.mass << " kg" << std::endl;
    std::cout << "LENGTH_UNIT: " << units.length << " m" << std::endl;
    std::cout << "VELOCITY_UNIT: " << units.velocity << " m/s" << std::endl;
    std::cout << "ENERGY_UNIT: " << units.energy << " J" << std::endl;
    std::cout << "TIME_UNIT: " << units.time << " s" << std::endl;
    std::cout << time << " " << end_time_nb << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    while (time<end_time_nb){
        leapfrog_kdk(time_step);
        time += time_step;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double total_energy = get_energy();
    double dx = pos[p1][0] - pos[p2][0];
    double dy = pos[p1][1] - pos[p2][1];
    double dz = pos[p1][2] - pos[p2][2];
    double dr2 = dx*dx + dy*dy + dz*dz;
    dr2 = length_nb_to_si(sqrt(dr2));

    std::cout << " r [au]: " << dr2 / 1.5e11 << std::endl;
    std::cout << "dE: " << (total_energy - initial_energy) / initial_energy << std::endl;
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " ms" << std::endl;

    return 0;
}