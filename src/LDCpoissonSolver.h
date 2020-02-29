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

private:

    // Methods
    void Build2DLaplace();
    void ComputeFactor();
    
    // Linear solver variables
    int nNodes, Nx, Ny, nCoeffs;
    int rank;

    double* A = nullptr;
    double coeff[3];
};