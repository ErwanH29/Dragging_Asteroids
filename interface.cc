/* 
Symplectic integrator with several choices of integrator.
Symplectic integrator_choices preserve energy error over 
secular times, however across an orbital period some 
error do arise. This error can be reduced by making the 
internal time-step an integer fraction of the orbital 
period, and further reduced by making it an integer 
fraction of the numerical orbital period, and even further 
reduced by making this a power of two integer fraction of 
orbital period (numerical round off).

NOTE: double is better for array and heavier objects, use
double otherwise.

Available integrator_choices:
    - Second-order verlet leapfrog method       22/05/2026
    - Fourth-order Yoshida integrator_choice    26/05/2026
    - Eighth-order Yoshida integrator_choice    26/05/2026
    - Wisdam-Holman + Bulirsch-Stoer
*/


#include <chrono>
#include <cmath>
#include <vector>

#include "nbody_system.h"
#include "orbital_elements.h"
#include "units.h"
#include "vec_3d.h"


int main(){
    double M_si = 2.0e30;  // Mass of the sun
    double r_si = 1.5e11;  // 1 AU
    int integrator_choice = 3;
    if (integrator_choice!=1 && integrator_choice!=2 && integrator_choice!=3){
        std::cout << "Invalid integrator choice. Current choices are 1, 2, and 3." << std::endl;
        return -1;
    }

    NBodySystem system;
    system.set_eps2(0.0);
    system.set_mass_threshold(pow(2.0, -15.0));

    system.new_particle(
        1.0, M_si, 
        0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0
    );
    system.new_particle(
        1.0, M_si, 
        -r_si, 0.0, 0.0, 
        0.0, sqrt(G_si * M_si / r_si), 0.0
    );
    for (int i=0; i<100; i++){
        double r = (0.8 + i/200.0) * r_si;
        system.new_particle(
            1.0, pow(2.0, -16.0), 
            -r, 0.0, 0.0, 
            0.0, sqrt(G_si * M_si / r), 0.0
        );
    }

    system.set_nb_units();
    system.convert_si_to_nb();

    system.update_gravity();
    system.move_to_center();
    const double initial_energy = system.get_energy();


if (1){
    std::cout << "MASS_UNIT: " << units.mass << " kg" << std::endl;
    std::cout << "LENGTH_UNIT: " << units.length << " m" << std::endl;
    std::cout << "VELOCITY_UNIT: " << units.velocity << " m/s" << std::endl;
    std::cout << "ENERGY_UNIT: " << units.energy << " J" << std::endl;
    std::cout << "TIME_UNIT: " << units.time << " s" << std::endl;
    std::cout << "==============" << std::endl;
}

    const double initial_ecc = get_ecc(
        system.mass[0], system.mass[1], 
        system.pos[0], system.pos[1], 
        system.vel[0], system.vel[1]
    );
    const double initial_sma = get_sma(
        system.mass[0], system.mass[1], 
        system.pos[0], system.pos[1], 
        system.vel[0], system.vel[1]
    );

    std::cout << "Initial ecc= " << initial_ecc << std::endl;
    std::cout << "Initial sma= " << initial_sma * units.length / 1.5e11 << " au" << std::endl;

    double orbital_period = get_orbital_period(
        system.mass[0], system.mass[1], initial_sma
    );


    int counter = 0;
    double time = 0.0;
    const double end_time_si = pow(10.0, 3.0) * 3600.0 * 24.0 * 365.25;
    const double end_time_nb = units.time_si_to_nb(end_time_si);
    const double time_step = orbital_period / pow(2.0, 5.0);
    const int n_steps = end_time_nb / time_step;
    std::cout << "Number of steps: " << n_steps << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    while (time<end_time_nb){
        if (integrator_choice==1){
            system.leapfrog_kdk(time_step);
        } else if (integrator_choice==2){
            system.yoshida_fourth_order(time_step);
        } else if (integrator_choice==3){
            system.yoshida_eighth_order(time_step);
        }
        time += time_step;
        counter++;
        if (counter%100000 == 0){
            std::cout << "\rProgress: " << (time/end_time_nb)*100 << "%" << std::flush;
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    system.update_gravity();


if (1){
    system.move_to_center();
    double total_energy = system.get_energy();
    double dx = system.pos[0][0] - system.pos[1][0];
    double dy = system.pos[0][1] - system.pos[1][1];
    double dz = system.pos[0][2] - system.pos[1][2];
    double dr2 = dx*dx + dy*dy + dz*dz;
    dr2 = units.length_nb_to_si(sqrt(dr2));

    std::cout << " r [au]: " << dr2 / 1.5e11 << std::endl;
    std::cout << "dE: " << (total_energy - initial_energy) / initial_energy << std::endl;
}


    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " ms" << std::endl;

    const double final_ecc = get_ecc(
        system.mass[0], system.mass[1], 
        system.pos[0], system.pos[1], 
        system.vel[0], system.vel[1]
    );
    const double final_sma = get_sma(
        system.mass[0], system.mass[1], 
        system.pos[0], system.pos[1], 
        system.vel[0], system.vel[1]
    );
    std::cout << "Final ecc= " << final_ecc << std::endl;
    std::cout << "Final sma= " << final_sma * units.length / 1.5e11 << " au" << std::endl;
    for (int i=0; i<system.size(); i++){
        std::cout << "Particle " << i << ": x-pos = " << system.pos[i][0] << " m" << std::endl;
    }
    return 0;
}