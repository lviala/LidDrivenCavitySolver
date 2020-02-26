#include <iostream>
#include <iterator>
#include <mpi.h>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "LDCprogram_options.h"
#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{   
    bool b_help = false; // Help menu request boolean
    int gridSize[2];
    double domainSize[2], timeStep, finalTime, reynoldsNumber;
    
    MPI_Init(&argc, &argv);

    ///////////////////////////////////////////////////////////////
    // MPI TEST
    int rank, size, retval_rank, retval_size;

    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank); // zero-based 
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM) {
        std::cout << "Invalid communicator" << std::endl;
    return 1; }

    cout << "I am process " << rank + 1 << " of " << size << endl;
    ///////////////////////////////////////////////////////////////

    // Read and parse program options from command line
    po::variables_map vm;
    b_help = LDCprogram_options(argc , argv, vm);

    // If help is displayed, return 0 and exit
    if (b_help) {
        cout << "Program complete with exit code 0" << endl;
        return 0;
    }

    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();
    
    // Assign problem parameters to solver
    LDCset_solver(solver , vm);
    
    // Initialize solver
    solver->Initialise();

    // Run the solver
    solver->Integrate();

    MPI_Finalize();

	return 0;
}