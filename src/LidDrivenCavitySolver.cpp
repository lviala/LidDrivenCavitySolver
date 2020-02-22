#include <iostream>
#include <iterator>
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
    solver -> getDomainSize(domainSize);
    cout << "Domain x-direction size: " << domainSize[0] << " meters" << endl;
    cout << "Domain y-direction size: " << domainSize[1] << " meters" << endl;

    solver -> getGridSize(gridSize);
    cout << "Grid x-direction size: " << gridSize[0] << " nodes" << endl;
    cout << "Grid y-direction size: " << gridSize[1] << " nodes" << endl;

    solver -> getTimeStep(timeStep);
    cout << "Time step: " << timeStep << " seconds" << endl;

    solver -> getFinalTime(finalTime);
    cout << "Final time: " << finalTime << " seconds" << endl;

    solver -> getReynoldsNumber(reynoldsNumber);
    cout << "Reynolds Number: " << reynoldsNumber << endl;
    

    // Initialize solver
    solver->Initialise();

    // Run the solver
    solver->Integrate();

	return 0;
}