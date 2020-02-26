#include "LidDrivenCavity.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

bool LDCprogram_options(int argc, char** argv, po::variables_map& vm, int& rank) {

    po::options_description desc("Allowed options");
    // Specify allowed options
    desc.add_options()
        ("help" , "Produce help message.")
        ("Lx" , po::value<double>() -> default_value(1.0)    , "Length of the domain in the x-direction.")
        ("Ly" , po::value<double>() -> default_value(1.0)    , "Length of the domain in the y-direction.")
        ("Nx" , po::value<int>() -> default_value(161)       , "Number of grid points in x-direction.")
        ("Ny" , po::value<int>() -> default_value(161)       , "Number of grid points in y-direction.")
        ("Px" , po::value<int>() -> default_value(1)         , "Number of partitions in x-direction.")
        ("Py" , po::value<int>() -> default_value(1)         , "Number of partitions in y-direction.")
        ("dt" , po::value<double>() -> default_value(0.01)   , "Time step size.")
        ("T"  , po::value<double>() -> default_value(10.0)   , "Final time.")
        ("Re" , po::value<double>() -> default_value(1000.0) , "Reynolds number.")
    ;

    // Parse command line arguments and store specified values
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    // If the user gives the --help argument, print the help and quit.
    if (vm.count("help")) {
        if (rank == 0){
            cout << desc << "\n";
        }
        return true;
    }
    else // Return false, and proceed with rest of program
    {
        return false;
    }
    

}

void LDCsetVar(po::variables_map& vm, 
                    int* gridSize, int* partitionSize, double* domainSize, 
                    double& timeStep, double& finalTime, double& reynoldsNumber) {

    // Note: working in column-major format: y-direction to be 0-direction
    gridSize[1] = vm["Nx"].as<int>();
    gridSize[0] = vm["Ny"].as<int>();

    partitionSize[1] = vm["Px"].as<int>();
    partitionSize[0] = vm["Py"].as<int>();

    domainSize[1] = vm["Lx"].as<double>();
    domainSize[0] = vm["Ly"].as<double>();

    timeStep = vm["dt"].as<double>();
    finalTime = vm["T"].as<double>();
    reynoldsNumber = vm["Re"].as<double>();

}