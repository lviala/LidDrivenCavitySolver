#pragma once

#include <string>
#include <cstring>

#include "LDCpoissonSolver.h"

using namespace std;

class LDCpoissonSolver_Banded : public LDCpoissonSolver
{
public:
    // Constructor
    LDCpoissonSolver_Banded(int rank):LDCpoissonSolver(rank){}
    ~LDCpoissonSolver_Banded();

    virtual void Initialize(int& Nx, int& Ny, double* coeff);
    virtual void SolvePoisson(double* v,double* s);

    // IO Methods
    virtual void PrintCoeffMat(int rank);
    void PrintRHS(int rank);

protected:

    // Methods
    void Build2DLaplace();
    void Factor2DLaplace();
    void BuildRHS(double* v, double* s);

};