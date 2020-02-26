#pragma once

#include <string>
using namespace std;

class LidDrivenCavity
{
public:
    // Constructors
    LidDrivenCavity(int rank, int* coords, int* gridSize, double dt, double T, double Re);
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
    
    // Methods
    void Initialise();
    void Integrate();

private:

    // Contains coordinate of subdomain on cartesian grid
    int coords[2];
    int rank;

    double* v = nullptr;
    double* s = nullptr;

    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Re;
};
