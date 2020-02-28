#pragma once

#include <string>
#include <cstring>
#include <mpi.h>
using namespace std;

class LidDrivenCavity
{
public:
    // Constructors
    LidDrivenCavity(MPI_Comm MPIcomm, int rank, int* rankShift, int* coords, int* gridSize, double dt, double dx, double dy, double T, double Re);
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
    
    // Solver Methods
    void Initialise();
    void UpdateGlobalBcs();
    void Integrate();

    // Interface Management Methods
    void InterfaceBroadcast(double* field);
    void InterfaceGather(double* field);
    
    // IO Methods
    void PrintArray(const char* varStr, int rank);
    void LDCStatus(int rank);

private:

    //Methods
    void InterfaceSend(int& count, double* field, double* buff, int disp, int& dest, int& tag, MPI_Comm MPIcomm);
    void InterfaceRecv(int& count, double* field, double* buff, int disp, int& source, int& tag, MPI_Comm MPIcomm);

    MPI_Comm MPIcomm;
    // Contains coordinate of subdomain on cartesian grid
    int coords[2];
    // rankShift[0] = rank of process below - rankShift[1] = rank of process above
    // rankShift[2] = rank of process left - rankShift[3] = rank of process right 
    int rankShift[4];
    int rank;

    // Streamfunction and Vorticity field arrays
    double* v = nullptr;
    double* s = nullptr;
    
    // Buffer arrays for MPI communications
    double* bufNx = nullptr;
    double* bufNy = nullptr;

    // Problem parameters
    double dt;
    double dx;
    double dy;
    double T;
    int    Nx;
    int    Ny;
    double Re;
    double U = 1.0;
};
