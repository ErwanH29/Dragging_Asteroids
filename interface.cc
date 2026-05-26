/* 
Symplectic integrator_choice with 1st, 2nd, 3rd and 4th
order schemes. Symplectic integrator_choices preserve
energy error over secular times, however across
an orbital period some error does arise. 

This error can be reduced by making the internal 
time-step an integer fraction of the orbital period. 
This error can be reduced further by making it an 
integer fraction of the numerical orbital period, 
and further reduced by making this a power of two 
integer fraction of orbital period (numerical round off).

NOTE: double is better for array and heavier objects, use
double otherwise.

Available integrator_choices:
    - Second-order verlet leapfrog method       22/05/2026
    - Fourth-order Yoshida integrator_choice
    - Eighth-order Yoshida integrator_choice
    - Wisdam-Holman + Bulirsch-Stoer
*/

#include <chrono>
#include <cmath>
#include <vector>

#include "orbital_elements.h"
#include "units.h"
#include "vec_3d.h"

#define _TINY_ pow(2.0, -52.0)



template <typename T>
using dyn_arr = std::vector<T>;
static dyn_arr<double> mass, pot, radius;
static dyn_arr<int> type;
static dyn_arr<vec> acc, pos, vel;

constexpr int NDIM = 3;  // Known at compile time
static int integrator_choice = 2;
static double eps2 = 0.0;  // only interface.cc can use it, but global
static double MASS_THRESHOLD = pow(2.0, -10.0);



void set_nb_units(){
    double mass_units = 0.0;
    double length_units = 0.0;
    int nparticles = radius.size();
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
};

void convert_si_attr_to_nb(){
    int nparticles = radius.size();
    for (int i=0; i<nparticles; i++){
        mass[i] = mass_si_to_nb(mass[i]);
        radius[i] = length_si_to_nb(radius[i]);
        pos[i] = pos[i] * (1.0 / units.length);
        vel[i] = vel[i] * (1.0 / units.velocity);
    }
}

double get_energy(){
    double ke = 0.0;
    double pe = 0.0;
    for (int i=0; i<radius.size(); i++){
        double mass_i = mass[i];
        if (mass_i <= _TINY_) continue;
        double vel2 = vel[i] * vel[i];  // Using overloaded * operator
        ke += 0.5 * mass_i * vel2;
        pe += 0.5 * mass_i * pot[i];
    }
    return ke + pe;
};


int new_particle( 
    double _radius,
    double _mass, 
    double x, 
    double y, 
    double z, 
    double vx, 
    double vy, 
    double vz
){
    int id = radius.size();
    static const vec null_acc = vec(0.0, 0.0, 0.0);

    mass.push_back(_mass);
    radius.push_back(_radius);
    pos.push_back(vec(x, y, z));
    vel.push_back(vec(vx, vy, vz));
    acc.push_back(null_acc);
    pot.push_back(0.0);
    if (_mass >= MASS_THRESHOLD){
        type.push_back(1);
    } else {
        type.push_back(0);
    }
    return id;
};


int get_acc_pot_coll(){
    int nparticles = radius.size();
    static vec zero_vec = vec(0.0, 0.0, 0.0);  // Initialised at first call, and remembered throughout
    for (int i=0; i<nparticles; i++){
        acc[i] = zero_vec;
        pot[i] = 0.0;
    }

    for (int i=0; i<nparticles; i++){
        const vec& pos_i = pos[i];
        const vec& vel_i = vel[i];
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
            vec vij = vel_i - vel[j];
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
};


int kick(double time_step){
    int nparticles = radius.size();
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            vel[i][k] += time_step * acc[i][k];
        }
    }
    return 0;
}


int main(){
    double M_si = 2.0e30;  // Mass of the sun
    double r_si = 1.5e11;  // 1 AU
    int (*integrator)(double) = nullptr;
    if (integrator_choice == 1){
        integrator = &leapfrog_kdk;
    }
    else if (integrator_choice == 2){
        integrator = &yoshida_fourth_order;
    }
    else{
        std::cout << "Invalid integrator choice. Current choices are 1, 2, and 3." << std::endl;
        return -1;
    }


    int p1 = new_particle(
        1.0, M_si, 
        0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0
    );
    int p2 = new_particle(
        1.0, 0.2 * M_si, 
        -r_si, 0.0, 0.0, 
        0.0, sqrt(G_si * M_si / r_si), 0.0
    );

    for (int i=0; i<10; i++){
        double r = (0.9 + i/2000.0) * r_si;
        int p = new_particle(
            1.0, pow(2.0, -16.0), 
            -r, 0.0, 0.0, 
            0.0, sqrt(G_si * M_si / (r)), 0.0
        );
    }

    set_nb_units();
    convert_si_attr_to_nb();

    get_acc_pot_coll();
    const double initial_energy = get_energy();

    std::cout << "MASS_UNIT: " << units.mass << " kg" << std::endl;
    std::cout << "LENGTH_UNIT: " << units.length << " m" << std::endl;
    std::cout << "VELOCITY_UNIT: " << units.velocity << " m/s" << std::endl;
    std::cout << "ENERGY_UNIT: " << units.energy << " J" << std::endl;
    std::cout << "TIME_UNIT: " << units.time << " s" << std::endl;
    std::cout << "==============" << std::endl;

    const double initial_ecc = get_ecc(
        mass[p1], mass[p2], 
        pos[p1], pos[p2], 
        vel[p1], vel[p2]
    );
    const double initial_sma = get_sma(
        mass[p1], mass[p2], 
        pos[p1], pos[p2], 
        vel[p1], vel[p2]
    );
    std::cout << "Initial ecc= " << initial_ecc << std::endl;
    std::cout << "Initial sma= " << initial_sma * units.length / 1.5e11 << " au" << std::endl;

    double orbital_period = get_orbital_period(
        mass[p1], mass[p2], initial_sma
    );


    double time = 0.0;
    const double end_time_si = 1000.0 * 3600.0 * 24.0 * 365.25;
    const double end_time_nb = time_si_to_nb(end_time_si);
    const double time_step = orbital_period / pow(2.0, 6.0);
    std::cout << "Number of steps: " << end_time_nb / time_step << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    while (time<end_time_nb){
        integrator(time_step);
        time += time_step;
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    get_acc_pot_coll();
    double total_energy = get_energy();
    double dx = pos[p1][0] - pos[p2][0];
    double dy = pos[p1][1] - pos[p2][1];
    double dz = pos[p1][2] - pos[p2][2];
    double dr2 = dx*dx + dy*dy + dz*dz;
    dr2 = length_nb_to_si(sqrt(dr2));

    std::cout << " r [au]: " << dr2 / 1.5e11 << std::endl;
    std::cout << "dE: " << (total_energy - initial_energy) / initial_energy << std::endl;
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " ms" << std::endl;

    const double final_ecc = get_ecc(
        mass[p1], mass[p2], 
        pos[p1], pos[p2], 
        vel[p1], vel[p2]
    );
    const double final_sma = get_sma(
        mass[p1], mass[p2], 
        pos[p1], pos[p2], 
        vel[p1], vel[p2]
    );
    std::cout << "Final ecc= " << final_ecc << std::endl;
    std::cout << "Final sma= " << final_sma * units.length / 1.5e11 << " au" << std::endl;
    return 0;
}