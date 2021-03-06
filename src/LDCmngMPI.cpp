#include "LidDrivenCavity.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
using namespace std;

namespace mngMPI{

    bool validateNP(int& np, int* partitionSize) {
        // Check if number of processes is compatible with number
        // of domain partitions

        if (np != partitionSize[0]*partitionSize[1]){
            return false;
        }
        else{
            return true;
        }

    }

    void splitGrid(int* gridSize, int* partitionSize, int* coords, int* subGridSize, double& dx, double& dy, double* subPos){
        // Splits grip into specified number of subdomains, distributing nodes as equally as possible
        int subNx_quotient, subNx_remainder, subNy_quotient, subNy_remainder;

        subNx_quotient = gridSize[1] / partitionSize[1];
        subNx_remainder = gridSize[1] % partitionSize[1];

        subNy_quotient = gridSize[0] / partitionSize[0];
        subNy_remainder = gridSize[0] % partitionSize[0];

        if (coords[1] <  subNx_remainder ){
            subGridSize [1] = subNx_quotient + 1;
            
            subPos[1] = dx * subGridSize [1] * coords[1];
        }
        else{
            subGridSize [1] = subNx_quotient;
            
            subPos[1] = dx * ((subGridSize [1] + 1) * subNx_remainder +
                            subGridSize [1] * (coords[1] - subNx_remainder));
        }

        if (coords[0] <  subNy_remainder ){
            subGridSize [0] = subNy_quotient + 1;
            
            subPos[0] = dy * subGridSize [0] * coords[0];
        }
        else{
            subGridSize [0] = subNy_quotient;

            subPos[0] = dy * ((subGridSize [0] + 1) * subNy_remainder +
                            subGridSize [0]  * (coords[0] - subNy_remainder));
        }
    }
}