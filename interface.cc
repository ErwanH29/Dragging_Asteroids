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

#include <cmath>
#include <vector>
#include "vec_3d.h"

#define _TINY_ pow(2.0, -52.0)

template <typename T>
using dyn_arr = std::vector<T>;

const double G = 6.67430e-11;
const int NDIM = 3;
static const double eps2 = 0.0;
static dyn_arr<double> mass, pot, radius;
static dyn_arr<vec> acc, pos, vel;


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

            double rij[NDIM], vij[NDIM];
            vec pos_j = pos[j];
            vec vel_j = vel[j];
            for (int k = 0; k<NDIM; k++){
                rij[k] = pos_i[k] - pos_j[k];
                vij[k] = vel_i[k] - vel_j[k];
            }

            double dx = pos_i[0] - pos_j[0];
            double dy = pos_i[1] - pos_j[1];
            double dz = pos_i[2] - pos_j[2];
            double dr2 = dx*dx + dy*dy + dz*dz;
            double coll_radius = rad_i + radius[j];
            if (dr2 <= coll_radius * coll_radius){
                // Handle collision (not yet implemented)
                continue;
            }

            dr2 += eps2;  // Softening to avoid singularities
            double dr = sqrt(dr2);
            double dr3 = dr * dr2;

            //Compute acc. and pot.
            double da[NDIM];
            for (int k=0; k<NDIM; k++){
                da[k] = -G * rij[k] / dr3;
            }
            if (mass_j>_TINY_){
                pot[i] -= G * mass_j / dr;
                for (int k=0; k<NDIM; k++){
                    acc[i][k] += mass_j * da[k];
                }
            }
            if (mass_i>_TINY_){
                pot[j] -= G * mass_i / dr;
                for (int k=0; k<NDIM; k++){
                    acc[j][k] -= mass_i * da[k];  // drji = -drij
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
    double M = 2e30;                // Mass of the sun
    double r = 1.5e11;              // 1 AU
    double vorb = sqrt(G * M / r);  // Orbital velocity
    int p1 = new_particle(
        1.0, M, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    );
    int p2 = new_particle(
        1.0, M, -r, 0.0, 0.0, 0.0, -0.99 * vorb, 0.0
    );

    double time = 0.0;
    double end_time = 100.0 * 3600 * 24 * 365.25;  // 100 years
    const double time_step = end_time / pow(2.0, 15.0);  // 32k steps per orbit

    get_acc_pot_coll();
    const double initial_energy = get_energy();
    while (time<end_time){
        leapfrog_kdk(time_step);

        double total_energy = get_energy();
        std::cout << "dE: " << (total_energy - initial_energy) / initial_energy;

        double dx = pos[p1][0] - pos[p2][0];
        double dy = pos[p1][1] - pos[p2][1];
        double dz = pos[p1][2] - pos[p2][2];
        double dr2 = dx*dx + dy*dy + dz*dz;
        std::cout << "r: " << sqrt(dr2) / 1.5e11 << std::endl;


        time += time_step;
    }

    return 0;
}