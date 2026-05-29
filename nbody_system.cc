#include <cmath>
#include "nbody_system.h"
#include "units.h"
#include "vec_3d.h"

NBodySystem::NBodySystem(){}

int NBodySystem::size()const{
    return radius.size();
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
            length_units = std::max(length_units, dr);
        }
    }
    units.set_units(mass_units, length_units);
}

void NBodySystem::convert_si_to_nb(){
    int nparticles = size();
    for (int i=0; i<nparticles; i++){
        mass[i] = units.mass_si_to_nb(mass[i]);
        radius[i] = units.length_si_to_nb(radius[i]);
        pos[i] = pos[i] * (1.0 / units.length);
        vel[i] = vel[i] * (1.0 / units.velocity);
    }
}

// Particle manipulation
void NBodySystem::new_particle(
    double _radius,
    double _mass, 
    double x, 
    double y, 
    double z, 
    double vx, 
    double vy, 
    double vz
){
    static const vec null_acc = vec(0.0, 0.0, 0.0);
    if (_mass <= _TINY_){
        _mass = 0.0;
    }
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
}

void NBodySystem::remove_particle(int i){
    radius.erase(radius.begin() + i);
    mass.erase(mass.begin() + i);
    pos.erase(pos.begin() + i);
    vel.erase(vel.begin() + i);
    acc.erase(acc.begin() + i);
    pot.erase(pot.begin() + i);
    type.erase(type.begin() + i);
}



// Setters
void NBodySystem::set_collision_detection(bool _collision_detection){
    collision_detection = _collision_detection;
}

void NBodySystem::set_eps2(double _eps2){
    eps2 = _eps2;
}

void NBodySystem::set_eta(double _eta){
    eta = _eta;
}

void NBodySystem::set_integrator_choice(int _integrator_choice){
    integrator_choice = _integrator_choice;
}

void NBodySystem::set_mass_threshold(double _mass_threshold){
    mass_threshold = _mass_threshold;
}


// Getters
bool NBodySystem::get_collision_detection() const {
    return collision_detection;
}

double NBodySystem::get_eps2() const {  //When reading, use const reference to avoid unnecessary copying
    return eps2;
}

double NBodySystem::get_eta() const {
    return eta;
}

double NBodySystem::get_mass_threshold() const {
    return mass_threshold;
}

int NBodySystem::get_integrator_choice() const {
    return integrator_choice;
}

int NBodySystem::get_ncoll() const {
    return coll_counter;
}


// Helper functions
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

void NBodySystem::move_to_center(){
    static vec zero_vec = vec(0.0, 0.0, 0.0);
    vec center_of_mass_pos = zero_vec;
    vec center_of_mass_vel = zero_vec;
    
    double total_mass = 0.0;

    int nparticles = size();
    for (int i=0; i<nparticles; i++){
        const double mass_i = mass[i];
        if (mass[i] <= _TINY_) continue;
        total_mass += mass[i];
        center_of_mass_pos += pos[i] * mass_i;
        center_of_mass_vel += vel[i] * mass_i;
    }

    center_of_mass_pos = center_of_mass_pos * (1.0 / total_mass);
    center_of_mass_vel = center_of_mass_vel * (1.0 / total_mass);

    for (int i=0; i<nparticles; i++){
        pos[i] -= center_of_mass_pos;
        vel[i] -= center_of_mass_vel;
    }
}


// Handle collisions
void NBodySystem::resolve_collision(int i, int j){
    double new_radius = std::cbrt(  // Assume constant density
        radius[i]*radius[i]*radius[i] +
        radius[j]*radius[j]*radius[j]
    );
    double new_mass = mass[i] + mass[j];
    vec new_vel = (mass[i]*vel[i] + mass[j]*vel[j]) / new_mass;
    vec new_pos = (mass[i]*pos[i] + mass[j]*pos[j]) / new_mass;

    NBodySystem::new_particle(
        new_radius,
        new_mass,
        new_pos[0],
        new_pos[1],
        new_pos[2],
        new_vel[0],
        new_vel[1],
        new_vel[2]
    );

    if (i>j){
        NBodySystem::remove_particle(i);
        NBodySystem::remove_particle(j);
    }
    else{
        NBodySystem::remove_particle(j);
        NBodySystem::remove_particle(i);
    }

    coll_counter++;
}

bool NBodySystem::detect_collision(){
    const int nparticles = size();

    for (int i = 0; i < nparticles; i++){
        for (int j = i + 1; j < nparticles; j++){

            // Only check collisions involving at least one massive particle
            if (mass[i] < mass_threshold && mass[j] < mass_threshold){
                continue;
            }

            const double dx = pos[i][0] - pos[j][0];
            const double dy = pos[i][1] - pos[j][1];
            const double dz = pos[i][2] - pos[j][2];

            const double dr2 = dx*dx + dy*dy + dz*dz;
            const double coll_rad = radius[i] + radius[j];

            if (dr2 <= coll_rad * coll_rad){
                resolve_collision(i, j);

                std::cout << "COLLISION " << coll_counter
                          << " DETECTED: Nparticles: "
                          << size() << std::endl;
                std::cout << "NEW MASS: " << mass.back() << std::endl;

                return true;
            }
        }
    }

    return false;
}



// Gravity calculations
void NBodySystem::interact_pair(int i, int j){
        const double dx = pos[i][0] - pos[j][0];
        const double dy = pos[i][1] - pos[j][1];
        const double dz = pos[i][2] - pos[j][2];

        double dr2 = dx*dx + dy*dy + dz*dz;

        // std::cout << "Mass[i]: " << mass[i] << std::endl;
        // std::cout << "Mass[j]: " << mass[j] << std::endl;

        dr2 += eps2;

        const double dr = std::sqrt(dr2);
        const double inv_dr = 1.0 / dr;
        const double inv_dr3 = inv_dr / dr2;

        const double fac = -G * inv_dr3;

        const double mi = mass[i];
        const double mj = mass[j];

        if (mj > _TINY_){
            pot[i] -= G * mj * inv_dr;

            const double mj_fac = mj * fac;
            acc[i][0] += mj_fac * dx;
            acc[i][1] += mj_fac * dy;
            acc[i][2] += mj_fac * dz;
        }

        if (mi > _TINY_){
            pot[j] -= G * mi * inv_dr;

            const double mi_fac = mi * fac;
            acc[j][0] -= mi_fac * dx;
            acc[j][1] -= mi_fac * dy;
            acc[j][2] -= mi_fac * dz;
        }
}

int NBodySystem::update_gravity(){
    if (collision_detection){
        while (detect_collision()){
            //Resolve collision in step until no more collisions
        }
    }
    const int nparticles = size();

    static const vec zero_vec = vec(0.0, 0.0, 0.0);
    for (int i = 0; i < nparticles; i++){
        acc[i] = zero_vec;
        pot[i] = 0.0;
    }

    std::vector<int> massive;
    std::vector<int> massless;
    
    massive.reserve(nparticles);
    massless.reserve(nparticles);

    for (int i = 0; i < nparticles; i++){
        if (type[i] == 1){
            massive.push_back(i);
        } else {
            massless.push_back(i);
        }
    }

    // Massive-massive interactions
    for (std::size_t a = 0; a < massive.size(); a++){
        for (std::size_t b = a + 1; b < massive.size(); b++){
            interact_pair(massive[a], massive[b]);
        }
    }

    // Massive-massless interactions
    for (int i : massive){
        for (int j : massless){
            interact_pair(i, j);
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
    // move_to_center();
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

    static const double c1 = 0.5 * w1;
    static const double c2 = 0.5 * (w0 + w1);
    static const double c3 = 0.5 * (w0 + w1);
    static const double c4 = w1 / 2.0;
    static const double d1 = w1;
    static const double d2 = w0;
    static const double d3 = w1;

    // move_to_center();
    drift(c1 * time_step);
    
    update_gravity();
    kick(d1 * time_step);
    drift(c2 * time_step);

    update_gravity();
    kick(d2 * time_step);
    drift(c3 * time_step);

    update_gravity();
    kick(d3 * time_step);
    drift(c4 * time_step);

    return 0;
}

int NBodySystem::yoshida_eighth_order(double time_step){
    // Coefficients from 1990PhLA..150..262Y
    static constexpr double w1 = -1.61582374150097;
    static constexpr double w2 = -2.44699182370524;
    static constexpr double w3 = -0.00716989419708120;
    static constexpr double w4 = 2.44002732616735;
    static constexpr double w5 = 0.157739928123617;
    static constexpr double w6 = 1.82020630970714;
    static constexpr double w7 = 1.04242620869991;
    static constexpr double w0 = 1.0 - 2.0 * (w1 + w2 + w3 + w4 + w5 + w6 + w7);

    static constexpr double w_seq[15] = {
        w7, w6, w5, w4, w3, w2, w1, w0, w1, w2, w3, w4, w5, w6, w7
    };

    // move_to_center();
    drift(0.5 * w_seq[0] * time_step);
    for (int i = 0; i < 15; i++) {

        update_gravity();
        kick(w_seq[i] * time_step);
        if (i < 14){
            drift(0.5 * (w_seq[i] + w_seq[i+1]) * time_step);
        }
        else {
            drift(0.5 * w_seq[i] * time_step);
        }
    }

    return 0;
}