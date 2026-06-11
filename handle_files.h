#ifndef HANDLE_FILES_CC
#define HANDLE_FILES_CC

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "handle_files.h"
#include "nbody_system.h"
#include "units.h"

using row_struct = std::vector<double>;  //Handle rows
using data_struct = std::vector<row_struct>;  //Handle data file


void record_snapshot(
    const NBodySystem& system, 
    double time_in_si,
    const std::string& filename
);
void record_collision(std::string filename);
data_struct read_file(const std::string filename);

#endif