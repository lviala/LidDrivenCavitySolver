#pragma once

#include <string>
using namespace std;

class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double Lx, double Ly);
    void SetGridSize(int Nx, int Ny);
    void SetTimeStep(double dt);
    void SetFinalTime(double T);
    void SetReynoldsNumber(double Re);

    void Initialise();
    void Integrate();

    // Add any other public functions

private:
    double* v = nullptr;
    double* s = nullptr;

    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
};
