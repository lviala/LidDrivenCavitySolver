#pragma once

#include <string>
#include <cstring>


using namespace std;

class LDCpoissonSolver
{
public:
    // Constructor
    LDCpoissonSolver(int rank);
    ~LDCpoissonSolver();

    void Initialize(int& Nx, int& Ny, double& dx, double& dy);

    // IO Methods
    void PrintTRIUArray(int rank);

private:

    // Methods
    void Build2DLaplace();
    void Factor2DLaplace();
    
    // Linear solver variables
    int nNodes, Nx, Ny, nCoeffs;
    int rank, info;

    double* A = nullptr;
    double coeff[3];
};