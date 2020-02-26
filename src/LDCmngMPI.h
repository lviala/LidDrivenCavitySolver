#pragma once

#include "LidDrivenCavity.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
using namespace std;
namespace mngMPI{
    
    bool validateNP(int& np, int* partitionSize);

    void splitGrid(int* gridSize, int* partitionSize, int* coords, int* subGridSize);

}
