#include "LidDrivenCavity.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

bool LDCprogram_options(int argc, char** argv, po::variables_map& vm) {

    po::options_description desc("Allowed options");
    // Specify allowed options
    desc.add_options()
        ("help" , "Produce help message.")
        ("Lx" , po::value<int>() -> default_value(1)         , "Length of the domain in the x-direction.")
        ("Ly" , po::value<int>() -> default_value(1)         , "Length of the domain in the y-direction.")
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
        cout << desc << "\n";
        return true;
    }
    else // Return false, and proceed with rest of program
    {
        return false;
    }
    

}

void LDCset_solver(LidDrivenCavity* solver , po::variables_map& vm) {

    solver -> SetDomainSize(vm["Lx"].as<int>() , vm["Ly"].as<int>() );
    solver -> SetGridSize(vm["Nx"].as<int>() , vm["Ny"].as<int>() );
    solver -> SetTimeStep(vm["dt"].as<double>() );
    solver -> SetFinalTime(vm["T"].as<double>() );
    solver -> SetReynoldsNumber(vm["Re"].as<double>() );

}