#include "LidDrivenCavity.h"

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
}

// Setters **************************************************
void LidDrivenCavity::SetDomainSize(double Lx, double Ly)
{
    this -> Lx = Lx;
    this -> Ly = Ly; 
}

void LidDrivenCavity::SetGridSize(int Nx, int Ny)
{
    this -> Nx = Nx;
    this -> Ny = Ny;
}

void LidDrivenCavity::SetTimeStep(double dt)
{
    this -> dt = dt;
}

void LidDrivenCavity::SetFinalTime(double T)
{
    this -> T = T;
}

void LidDrivenCavity::SetReynoldsNumber(double Re)
{
    this -> Re = Re;
}

// Getters **************************************************
void LidDrivenCavity::getDomainSize(double* domainSize)
{
    domainSize[0] = this -> Lx;
    domainSize[1] = this -> Ly;
}

void LidDrivenCavity::getGridSize(int* gridSize)
{
    gridSize[0] = this -> Nx;
    gridSize[1] = this -> Ny;
}

void LidDrivenCavity::getTimeStep(double& dt)
{
    dt = this -> dt;
}

void LidDrivenCavity::getFinalTime(double& T)
{
    T = this -> T;
}

void LidDrivenCavity::getReynoldsNumber(double& Re)
{
    Re = this -> Re;
}

// Solvers **************************************************
void LidDrivenCavity::Initialise()
{
}

void LidDrivenCavity::Integrate()
{
}