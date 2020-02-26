#pragma once

using namespace std;
#include <boost/program_options.hpp>

int LDCprogram_options(int argc, char** argv, po::variables_map& vm, int& rank);
void LDCsetVar(po::variables_map& vm, int* gridSize, int* partitionSize, double* domainSize, double& timeStep, double& finalTime, double& reynoldsNumber);