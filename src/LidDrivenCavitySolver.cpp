#include <iostream>
#include <iterator>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "LDCprogram_options.h"
#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    // Read and parse program options from command line
    po::options_description desc("Allowed options");
    LDCprogram_options(argc , argv, desc);

    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver here...
    // ...
    solver->Initialise();

    // Run the solver
    solver->Integrate();

	return 0;
}