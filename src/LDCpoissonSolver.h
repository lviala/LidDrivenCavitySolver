#pragma once

#include <string>
#include <cstring>


using namespace std;

class LDCpoissonSolver
// Abstract class implement common logic to solving linear system
// Solving algorithm specific methods are left as virtual
// Allows flexibility to implement different algorithms as needed
{
public:
    // Constructor
    LDCpoissonSolver(int rank);
    virtual ~LDCpoissonSolver() = default;

    // Solver Methods
    virtual void Initialize(int& Nx, int& Ny, double* coeff) = 0;
    virtual void SolvePoisson(double* v,double* s) = 0;

    // IO Methods
    void PrintRHS(int rank);
    virtual void PrintCoeffMat(int rank) = 0;
    
protected:

    // Solver Methods
    void BuildRHS(double* v, double* s);
    
    // Linear solver variables
    int nNodes, Nx, Ny, nCoeffs;
    int rank, info;

    double* A = nullptr;
    double* b = nullptr;
    double coeff[3];
};