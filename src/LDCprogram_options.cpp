#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

void LDCprogram_options(int argc, char** argv, po::options_description& desc) {

    // Specify allowed options
    desc.add_options()
        ("help" , "Produce help message.")
        ("Lx" , po::value<int>() , "Length of the domain in the x-direction.")
        ("Ly" , po::value<int>() , "Length of the domain in the y-direction.")
        ("Nx" , po::value<int>() , "Number of grid points in x-direction.")
        ("Ny" , po::value<int>() , "Number of grid points in y-direction.")
        ("Px" , po::value<int>() , "Number of partitions in x-direction.")
        ("Py" , po::value<int>() , "Number of partitions in y-direction.")
        ("dt" , po::value<double>() , "Time step size.")
        ("T" , po::value<double>() , "Final time.")
        ("Re" , po::value<double>() , "Reynolds number.");

    // Parse command line arguments and store specified values
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
        
    // If the user gives the --help argument, print the help and quit.
    if (vm.count("help")) {
        cout << desc << "\n";
    }

}