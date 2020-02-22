#include "LidDrivenCavity.h"

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
}

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

void LidDrivenCavity::Initialise()
{
}

void LidDrivenCavity::Integrate()
{
}