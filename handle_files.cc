#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "nbody_system.h"
#include "units.h"

using row_struct = std::vector<double>;  //Handle rows
using data_struct = std::vector<row_struct>;  //Handle data file


bool first_write = true;  // To control whether to write header in snapshot file

/*
Files are assumed to host:
    snapshot_time
    particle_id, particle_mass, particle_pos, particle_vels
*/
void record_snapshot(
    const NBodySystem& system, //Reference to avoid copying
    double time_in_si,
    const std::string& filename
){
    std::ofstream file;
    if (first_write){
        file.open(filename);
        first_write = false;
    } else {
        file.open(filename, std::ios::app);
    }

    if (!file.is_open()){
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    file << std::setprecision(10);
    file << "snapshot_time: " << time_in_si * s_to_yr << " yr\n";
    file << "particle_id, particle_mass, particle_pos_x, particle_pos_y, particle_pos_z, particle_vel_x, particle_vel_y, particle_vel_z\n";
    for (int i=0; i<system.size(); i++){
        file << i << ", " << system.mass[i] << ", " 
             << units.length_nb_to_si(system.pos[i][0]) << ", " 
             << units.length_nb_to_si(system.pos[i][1]) << ", " 
             << units.length_nb_to_si(system.pos[i][2]) << ", "
             << units.velocity_nb_to_si(system.vel[i][0])  << ", " 
             << units.velocity_nb_to_si(system.vel[i][1])  << ", " 
             << units.velocity_nb_to_si(system.vel[i][2])  << "\n";
    }

};


void record_collision(std::string filename){

};

data_struct read_file(const std::string filename){
    data_struct data;
    std::ifstream file(filename);
    if (!file.is_open()){
        std::cerr << "Error opening file: " << filename << std::endl;
        return data;
    }

    std::string line; // Temporarily store line from file
    while (std::getline(file, line)){  // Store each line in file within `line`
        row_struct row;
        std::stringstream ss(line);  // Allows to parse line
        std::string cell;  // Store CSV cell value
        while (std::getline(ss, cell, ',')){
            row.push_back(std::stod(cell));
        }
        data.push_back(row);
    }

    return data;
}
