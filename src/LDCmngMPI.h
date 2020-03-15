#pragma once

#include "LidDrivenCavity.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
using namespace std;
namespace mngMPI{
    // Build and validate parrallel subdomain grids
    bool validateNP(int& np, int* partitionSize);

    void splitGrid(int* gridSize, int* partitionSize, int* coords, int* subGridSize, double& dx, double& dy, double* subPos);

}
