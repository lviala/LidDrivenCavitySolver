#pragma once

#include <string>
#include <cstring>
#include <mpi.h>

#include "LDCpoissonSolver.h"
using namespace std;

class LidDrivenCavity
{

public:

    // Constructors
    LidDrivenCavity(MPI_Comm MPIcomm, int rank, int* rankShift, int* coords, int* gridSize, double dt, double dx, double dy, double* subPos, double T, double Re);
    ~LidDrivenCavity();

    // Setters
    void SetCoords(int* coords);
    void SetGridSize(int Nx, int Ny);
    void SetTimeStep(double dt);
    void SetFinalTime(double T);
    void SetReynoldsNumber(double Re);

    // Getters
    void getCoords(int* coords);
    void getGridSize(int* gridSize);
    void getTimeStep(double& dt);
    void getFinalTime(double& T);
    void getReynoldsNumber(double& Re);
    void getStreamFunction(double* s);
    void getVorticity(double* v);

    // Public solver Methods
    void Initialise();
    void Solve();

    // IO Methods
    void PrintArray(const char* varStr, int rank);
    void LDCStatus(int rank);
    void LDCPrintSolution2File(string filename, const double& wallTime);

private:
    // MEMBER CLASSES
    MPI_Comm MPIcomm;
    LDCpoissonSolver* poissonSolver;

    // Solver methods
    void UpdateGlobalBcs();
    void UpdateInteriorVorticity();
    void Integrate();

    // Finite Difference Operator Methods
    void FDLalplacianOperator(const double& alpha, double* x, double* y);
    void FDAdvectionOperator(const double& alpha, double* s_new, double* v_old, double* v_new);
    void FDCurlOperator(double alpha_x, double alpha_y, double* f, double* df_dx, double* df_dy);

    // Interface Management Methods
    void InterfaceBroadcast(double* field);
    void InterfaceGather(double* field);
    void InterfaceSend(int& count, double* field, double* buff, int disp, int& dest, int& tag, MPI_Comm MPIcomm);
    void InterfaceRecv(int& count, double* field, double* buff, int disp, int& source, int& tag, MPI_Comm MPIcomm);

    //MEMBER VARIABLES
    int coords[2];// Contains coordinate of subdomain on cartesian grid

    // rankShift[0] = rank of process below - rankShift[1] = rank of process above
    // rankShift[2] = rank of process left - rankShift[3] = rank of process right
    int rankShift[4];
    int rank;

    // Streamfunction and Vorticity field arrays
    double coeff[3];
    double* s = nullptr;
    double* v = nullptr;
    double* v_new = nullptr;
    double* test =nullptr;

    // Velocity Field Arrays
    double* velU = nullptr;
    double* velV = nullptr;

    // Buffer arrays for MPI communications
    double* bufNx = nullptr;
    double* bufNy = nullptr;

    // Problem parameters
    double dt;
    double dx;
    double dy;
    double subPos[2];
    double T;
    int    Nx;
    int    Ny;
    double Re;
    double U = 1.0;
};
