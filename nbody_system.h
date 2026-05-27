#ifndef NBODY_SYSTEM_H
#define NBODY_SYSTEM_H

#include <cmath>
#include <vector>
#include "vec_3d.h"

template <typename T>
using dyn_arr = std::vector<T>;

class NBodySystem{
    private:
        static constexpr int NDIM = 3;
        double eps2 = 0.0;
        double mass_threshold = pow(2.0, -16.0);
        int integrator_choice = 1; 

    public:

        dyn_arr<double> mass;
        dyn_arr<double> pot;
        dyn_arr<double> radius;
        dyn_arr<int> type;
        dyn_arr<vec> pos;
        dyn_arr<vec> vel;
        dyn_arr<vec> acc;

        NBodySystem();

        void set_eps2(double _eps2);
        double get_eps2();

        void set_mass_threshold(double _mass_threshold);
        double get_mass_threshold();

        void set_integrator_choice(int _integrator_choice);
        int get_integrator_choice();

        int size();

        int new_particle(
            double _radius,
            double _mass, 
            double x, 
            double y, 
            double z, 
            double vx, 
            double vy, 
            double vz
        );

        void set_nb_units();
        void convert_si_attr_to_nb();

        double get_energy();
        int update_gravity();
        int kick(double time_step);
        int drift(double time_step);

        int leapfrog_kdk(double time_step);
        int yoshida_fourth_order(double time_step);
        int yoshida_eighth_order(double time_step);

};

#endif