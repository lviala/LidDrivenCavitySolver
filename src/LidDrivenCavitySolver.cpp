#include <iostream>
#include <iterator>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "LDCprogram_options.h"
#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{   
    bool b_help = false;

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

    // Configure the solver here...
    // ...
    solver->Initialise();

    // Run the solver
    solver->Integrate();

	return 0;
}