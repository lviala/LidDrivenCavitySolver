#include <iostream>
#include <string>
#include <iterator>
#include <mpi.h>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "LDCprogram_options.h"
#include "LidDrivenCavity.h"
#include "LDCpoissonSolver_Packed.h"
#include "LDCmngMPI.h" // namespace mngMPI

int main(int argc, char **argv)
{

/////////////////////////////////////////////////////////////////////////////////////////////////
// INITIALIZE MPI, READ AND VALIDATE INPUTS
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // MPI related variables
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Domain related Variables
    bool b_help = false; // Help menu request boolean
    int gridSize[2], partitionSize[2];
    double domainSize[2], timeStep, xStep, yStep, finalTime, reynoldsNumber;

    // Read and parse program options from command line
    po::variables_map vm;
    b_help = LDCprogram_options(argc , argv, vm , rank);

    // If help is displayed, return 0 and exit
    if (b_help) {
        if(rank == 0) cout << "Program complete with exit code 0" << endl;
        MPI_Finalize();
        return 0;
    }

    // Pass argument values to main variables
    LDCsetVar(vm, gridSize, partitionSize, domainSize, timeStep, xStep, yStep, finalTime, reynoldsNumber );

    // If number of processes is not compatible with domain
    // partitions return 0 and exit;
    if (!mngMPI::validateNP(size, partitionSize)){

        if (rank == 0){
            cout << endl << "Please enter a valid combination of processes and domain partitions such that:" << endl
                    << "np = Px * Py" << endl << endl
                    << "Program complete with exit code 0" << endl;
        }

        MPI_Finalize();
        return 0;
    }

    // Validate time step restriction condition dt >= 0.25*dx*dy*Re
    // If invalid, return 0 and exit
    if (timeStep >= 0.25 * xStep * yStep * reynoldsNumber){
        if (rank == 0){
            cout << endl << "Invalid time step dt = " << timeStep << endl
                    << "Time step must satisfy the following condition:" << endl
                    << "dt >= 0.25 * dx * dy * Re" << endl
                    << "For the chosen values of dx, dy, Re, please enter a time step such that: " << endl
                    << "dt <=  " << 0.25 * xStep * yStep * reynoldsNumber << endl << endl
                    << "Program complete with exit code 0" << endl;
        }

        MPI_Finalize();
        return 0;
    }
    

/////////////////////////////////////////////////////////////////////////////////////////////////
// INITIALIZE SUBDOMAIN

    // Timing execution time
    double tStart, tEnd; 
    tStart = MPI_Wtime(); 

    // Initialize cartesian grid communicator
    MPI_Comm cartGrid;
    int periods[2] = {0, 0}, coords[2];
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, partitionSize, periods, reorder, &cartGrid);
    MPI_Cart_coords(cartGrid, rank, 2, coords);

    // Compute neighbor process ranks
    int rankShift[4] = {rank,rank,rank,rank};
    MPI_Cart_shift(cartGrid, 0, 1, &rankShift[0], &rankShift[1]);
    MPI_Cart_shift(cartGrid, 1, 1, &rankShift[2], &rankShift[3]);

    // Compute number of nodes assigned to subdomain
    int subGridSize[2];
    double subPos[2];
    mngMPI::splitGrid(gridSize, partitionSize, coords, subGridSize, xStep, yStep, subPos);

    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity(MPI_COMM_WORLD,
                                                rank, rankShift,
                                                coords, subGridSize,
                                                timeStep, xStep, yStep,
                                                subPos, finalTime, reynoldsNumber);

/////////////////////////////////////////////////////////////////////////////////////////////////
// SOLVE PROBLEM

    // Initialize solver
    solver->Initialise();

    // Run the solver
    solver->Solve();
    
    // End Time 
    tEnd = MPI_Wtime(); 

/////////////////////////////////////////////////////////////////////////////////////////////////
// OUTPUT SOLUTION AND EXIT

    // Output solution to file
    string resultsOutputPath;
    resultsOutputPath = string("results/LDCoutput") + 
                        string("_Lx_") + to_string(int(domainSize[1])) + string("_Ly_") + to_string(int(domainSize[0])) + 
                        string("_Nx_") + to_string(gridSize[1]) + string("_Ny_") + to_string(gridSize[0]) + 
                        string("_Px_") + to_string(partitionSize[1]) + string("_Py_") + to_string(partitionSize[0]) + 
                        string("_T_") + to_string(int(finalTime)) + 
                        string("_Re_") + to_string(int(reynoldsNumber)) + string(".csv");

    cout << "File output to: " << resultsOutputPath << endl;

    solver->LDCPrintSolution2File(resultsOutputPath, (tEnd - tStart));

    // Cleanup on program exit
    delete solver;
    MPI_Finalize();

	return 0;
}
