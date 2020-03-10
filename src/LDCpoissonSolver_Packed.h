#pragma once

#include <string>
#include <cstring>

#include "LDCpoissonSolver.h"

using namespace std;

class LDCpoissonSolver_Packed : public LDCpoissonSolver
{
public:
    // Constructor
    LDCpoissonSolver_Packed(int rank):LDCpoissonSolver(rank){}
    ~LDCpoissonSolver_Packed();

    virtual void Initialize(int& Nx, int& Ny, double* coeff);
    virtual void SolvePoisson(double* v,double* s);

    // IO Methods
    void PrintCoeffMat(int rank);
    void PrintRHS(int rank);

protected:

    // Methods
    void Build2DLaplace();
    void Factor2DLaplace();
    void BuildRHS(double* v, double* s);
    
    // // Linear solver variables
    // int nNodes, Nx, Ny, nCoeffs;
    // int rank, info;

    // double* A = nullptr;
    // double* b = nullptr;
    // double coeff[3];
};