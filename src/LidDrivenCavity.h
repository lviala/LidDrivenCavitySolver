#pragma once

#include <string>
#include <cstring>
using namespace std;

class LidDrivenCavity
{
public:
    // Constructors
    LidDrivenCavity(int rank, int* rankShift, int* coords, int* gridSize, double dt, double dx, double dy, double T, double Re);
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
    
    // Methods    
    void Initialise();
    void UpdateGlobalBcs();
    void Integrate();
    void PrintArray(const char* varStr, int rank);
    void LDCStatus(int rank);

private:

    // Contains coordinate of subdomain on cartesian grid
    int coords[2];
    // rankShift[0] = rank of process below - rankShift[1] = rank of process above
    // rankShift[2] = rank of process left - rankShift[3] = rank of process right 
    int rankShift[4];
    int rank;

    double* v = nullptr;
    double* s = nullptr;

    double dt;
    double dx;
    double dy;
    double T;
    int    Nx;
    int    Ny;
    double Re;
    double U = 1.0;
};
