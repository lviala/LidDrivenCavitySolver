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
}