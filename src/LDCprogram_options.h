#pragma once

using namespace std;
#include <boost/program_options.hpp>
#include "LidDrivenCavity.h"

int LDCprogram_options(int argc, char** argv, po::variables_map& vm);
void LDCset_solver(LidDrivenCavity* solver , po::variables_map& vm);