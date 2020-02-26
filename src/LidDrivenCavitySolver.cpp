#include <iostream>
#include <iterator>
#include <mpi.h>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "LDCprogram_options.h"
#include "LidDrivenCavity.h"
#include "LDCmngMPI.h"
//namespace mngMPI = LDCmngMPI;

int main(int argc, char **argv)
{   

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // MPI related variables
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Domain related Variables
    bool b_help = false; // Help menu request boolean
    int gridSize[2], partitionSize[2];
    double domainSize[2], timeStep, finalTime, reynoldsNumber;

    // Read and parse program options from command line
    po::variables_map vm;
    b_help = LDCprogram_options(argc , argv, vm);

    // If help is displayed, return 0 and exit
    if (b_help) {
        cout << "Program complete with exit code 0" << endl;
        MPI_Finalize();
        return 0;
    }

    // Pass argument values to main variables
    LDCsetVar(vm, gridSize, partitionSize, domainSize, timeStep, finalTime, reynoldsNumber );

    // If number of processes is not compatible with domain
    // domain partitions return 1 and exit;
    if (!mngMPI::validateNP(size, partitionSize)){
        
        if (rank == 0){
        cout << endl << "Please enter a valid combination of processes and domain partitions such that:" << endl
                << "np = Px * Py" << endl << endl
                << "Program complete with exit code 0" << endl;
        }

        MPI_Finalize();
        return 0;
    }

    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();
    
    // Initialize solver
    solver->Initialise();

    // Run the solver
    solver->Integrate();

    MPI_Finalize();

	return 0;
}