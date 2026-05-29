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

Still to do:
    - Parallelise computations
    - Resolving collisions in a more optimised fashion
    - Read AMUSE particle set
    - Output phase-space coordinates
    - Output collision information
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
    system.set_eta(pow(2.0, -7.0));
    system.set_eps2(0.0);
    system.set_mass_threshold(pow(2.0, -15.0));
    system.set_collision_detection(0);

    system.new_particle(  // Sun
        697500000, M_si, 
        0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0
    );
    system.new_particle( // Jupiter
        142984000 / 2.0, 
        1.89813e27, 
        5.4 * r_si, 0.0, 0.0, 
        0.0, sqrt(G_si * M_si / (5.4 * r_si)), 0.0
    );
    system.new_particle( // Saturn
        120536000 / 2.0, 
        5.683e26, 
        8.7 * r_si, 0.0, 0.0, 
        0.0, sqrt(G_si * M_si / (8.7 * r_si)), 0.0
    );
    system.new_particle( // Uranus
        51118000 / 2.0, 
        8.681e25, 
        15.0 * r_si, 0.0, 0.0, 
        0.0, sqrt(G_si * M_si / ( 15.0 * r_si)), 0.0
    );
    system.new_particle( // Neptune
        49528000 / 2.0, 
        1.024e26, 
        28.0 * r_si, 0.0, 0.0, 
        0.0, sqrt(G_si * M_si / (28.0 * r_si)), 0.0
    );
    for (int i=0; i<0; i++){
        double r = (5.35 + i/10000.0) * r_si; // 20 to 30 AU
        system.new_particle(
            1000000 / 2.0, 
            pow(2.0, -16.0), 
            r, 0.0, 0.0, 
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

    double min_sma = std::numeric_limits<double>::max();
    double target_mass;
    for (int i=1; i<5; i++){
        const double initial_sma = get_sma(
            system.mass[0], system.mass[i], 
            system.pos[0], system.pos[i], 
            system.vel[0], system.vel[i]
        );
        if (initial_sma < min_sma){
            min_sma = initial_sma;
            target_mass = system.mass[i];
        }
    }

if (1){
for (int i=1; i<5; i++){
    const double initial_ecc = get_ecc(
        system.mass[0], system.mass[i], 
        system.pos[0], system.pos[i], 
        system.vel[0], system.vel[i]
    );
    const double initial_sma = get_sma(
        system.mass[0], system.mass[i], 
        system.pos[0], system.pos[i], 
        system.vel[0], system.vel[i]
    );
    std::cout << "==============" << std::endl;
    std::cout << "Particle mass=" << system.mass[i] << ": " << std::endl;
    std::cout << "Initial ecc= " << initial_ecc << std::endl;
    std::cout << "Initial sma= " << initial_sma * units.length / 1.5e11 << " au" << std::endl;
}
}

    double orbital_period = get_orbital_period(
        system.mass[0], target_mass, min_sma
    );
    std::cout << "Orbital period: " << orbital_period * units.time / (3600.0 * 24.0 * 365.25) << " years" << std::endl;


    int counter = 0;
    double time = 0.0;
    const double end_time_si = pow(10.0, 5.0) * 3600.0 * 24.0 * 365.25;
    const double end_time_nb = units.time_si_to_nb(end_time_si);
    const double time_step = orbital_period * system.get_eta();
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
        if (counter%10000000 == 0){
            std::cout << "\rProgress: " << (time/end_time_nb)*100 << "%" << std::flush;
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    system.update_gravity();

    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " ms" << std::endl;


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
    std::cout << "Ncoll: " << system.get_ncoll() << std::endl;

    for (int i=1; i<5; i++){
        const double final_ecc = get_ecc(
            system.mass[0], system.mass[i], 
            system.pos[0], system.pos[i], 
            system.vel[0], system.vel[i]
        );
        const double final_sma = get_sma(
            system.mass[0], system.mass[i], 
            system.pos[0], system.pos[i], 
            system.vel[0], system.vel[i]
        );
        std::cout << "==============" << std::endl;
        std::cout << "Particle mass=" << system.mass[i] << ": " << std::endl;
        std::cout << "Final ecc= " << final_ecc << std::endl;
        std::cout << "Final sma= " << final_sma * units.length / 1.5e11 << " au" << std::endl;
    }
    for (int i=0; i<system.size(); i++){
        std::cout << "Particle " << i << ": x-pos = " << system.pos[i][0] << " m" << std::endl;
    }
}


    return 0;
}