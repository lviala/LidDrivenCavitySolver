#pragma once

#include <string>
#include <cstring>


using namespace std;

class LDCpoissonSolver
{
public:
    // Constructor
    LDCpoissonSolver(int rank);
    virtual ~LDCpoissonSolver() = default;

    virtual void Initialize(int& Nx, int& Ny, double* coeff) = 0;
    virtual void SolvePoisson(double* v,double* s) = 0;

    // IO Methods
    void PrintRHS(int rank);
    virtual void PrintCoeffMat(int rank) = 0;
    
protected:

    // Methods
    void BuildRHS(double* v, double* s);
    
    // Linear solver variables
    int nNodes, Nx, Ny, nCoeffs;
    int rank, info;

    double* A = nullptr;
    double* b = nullptr;
    double coeff[3];
};