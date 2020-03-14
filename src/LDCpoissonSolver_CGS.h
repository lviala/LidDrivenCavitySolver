#pragma once

#include <string>
#include <cstring>

#include "LDCpoissonSolver.h"

using namespace std;

class LDCpoissonSolver_CGS : public LDCpoissonSolver
{
public:
    // Constructor
    LDCpoissonSolver_CGS(int rank, int* rankShift, MPI_Comm MPIcomm)
                        :LDCpoissonSolver(rank)
                        {this -> MPIcomm = MPIcomm;
                         this -> rankShift[0] = rankShift[0];
                         this -> rankShift[1] = rankShift[1];
                         this -> rankShift[2] = rankShift[2];
                         this -> rankShift[3] = rankShift[3];}
    ~LDCpoissonSolver_CGS();

    virtual void Initialize(int& Nx, int& Ny, double* coeff);
    virtual void SolvePoisson(double* v,double* s);

    // IO Methods
    virtual void PrintCoeffMat(int rank);
    void PrintRHS(int rank);

protected:

    // Methods
    void Build2DLaplace();
    void BuildRHS(double* v, double* s);

    // Interface management method
    void InterfaceBroadcast(double* field);
    void InterfaceGather(double* field);
    void InterfaceSend(int& count, double* field, double* buff, int disp, int& dest, int& tag, MPI_Comm MPIcomm);
    void InterfaceRecv(int& count, double& alpha, double* field, double* buff, int disp, int& dest, int& tag, MPI_Comm MPIcomm);

    MPI_Comm MPIcomm;
    // Contains rank of neighbor processes
    int rankShift[4];

    double* Ap = nullptr;
    double* p = nullptr;
    double* r = nullptr;
    double* u = nullptr;

    double* matvecbuf_Nx = nullptr;
    double* matvecbuf_Ny = nullptr;

};