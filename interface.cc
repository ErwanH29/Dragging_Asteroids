// Symplectic integrator with 1st, 2nd, 3rd and 4th
// order schemes. Symplectic integrators preserve
// energy error over secular times, however across
// an orbital period some error does arise. 
//
// This error can be reduced by making the internal 
// time-step an integer fraction of the orbital period. 
// This error can be reduced further by making it an 
// integer fraction of the numerical orbital period, 
// and further reduced by making this a power of two 
// integer fraction of orbital period (numerical round off).

#include <cmath>
#include <vector>
#include "vec_3d.h"

#define _TINY_ pow(2.0, -52.0)

template <typename T>
using dyn_arr = std::vector<T>;

// Set-up initial arrays
static dyn_arr<int> identity;
static dyn_arr<double> mass, radius;
static dyn_arr<vec> pos, vel, acc, pot;

const int NDIM = 3;


int new_particle(
    int *id, 
    double _radius,
    double _mass, 
    double x, 
    double y, 
    double z, 
    double vx, 
    double vy, 
    double vz
){
    radius.push_back(_radius);
    mass.push_back(_mass);
    pos.push_back(vec(x, y, z));  //push_back takes one arg, so vec the three dimensions
    vel.push_back(vec(vx, vy, vz));
    acc.push_back(vec(0., 0., 0.));
    pot.push_back(vec(0., 0., 0.));

    return 0;
};


int get_acc_pot(){
    const int Nparticles = radius.size();

    for (int i = 0; i < Nparticles; i++){
        vec pos_i = pos[i];
        vec vel_i = vel[i];
        double rad_i = radius[i];
        for (int j = i+1; j < Nparticles; j++){
            double mass_j = mass[j];
            // Test particles don't affect one another
            if ((mass_j <= _TINY_) && (mass[i] <= _TINY_)){
                continue;
            }

            double rij[NDIM], vij[NDIM];
            vec pos_j = pos[j];
            vec vel_j = vel[j];
            for (int k = 0; k < NDIM; k++){
                rij[k] = pos_i[k] - pos_j[k];
                vij[k] = vel_i[k] - vel_j[k];
            }

            double dx = pos_i[0] - pos_j[0];
            double dy = pos_i[1] - pos_j[1];
            double dz = pos_i[2] - pos_j[2];

            double dr2 = dx*dx + dy*dy + dz*dz;
            double dr = sqrt(dr2);
            double dr3 = dr * dr2;
        }
    }
    return 0;
};

int first_order_leapfrog(){
    return 0;
}


int main(){
    new_particle(0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
    new_particle(0, 1.0, 1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0);

    double time = 0.0;
    double end_time = 100.0;
    double time_step = end_time / pow(2, 10);  // Example: 1024 time steps

    return 0;
}