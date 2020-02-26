#include "LidDrivenCavity.h"
#include <iostream>

LidDrivenCavity::LidDrivenCavity(int rank, int* rankShift, int* coords, int* gridSize, double dt, double T, double Re)
{   
    this -> rank = rank;
    this -> coords[0] = coords[0];
    this -> coords[1] = coords[1];
    this -> rankShift[0] = rankShift[0];
    this -> rankShift[1] = rankShift[1];
    this -> rankShift[2] = rankShift[2];
    this -> rankShift[3] = rankShift[3];

    this -> Nx = gridSize[1];
    this -> Ny = gridSize[0];
    
    this -> dt = dt;
    this -> T = T;
    this -> Re = Re;

    cout << "My rank is: " << this -> rank << endl; 
    cout << "My Coordinates are: (" << this -> coords[0] << " , " << this -> coords[1] << ")" << endl;
    cout << "My Neighbors are: Down=" << rankShift[0] << "  -- Up =" << rankShift[1] << "  -- Left =" << rankShift[2] << "  -- Right =" << rankShift[3] << endl;
    cout << "Nx=" << this -> Nx << " -- Ny=" << this -> Ny << endl;
}

LidDrivenCavity::~LidDrivenCavity()
{
}

//////////////////////////////////////////////////////////////
// SETTERS

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

void LidDrivenCavity::SetCoords(int* coords)
{
    this -> coords[0] = coords[0];
    this -> coords[1] = coords[1];
}

//////////////////////////////////////////////////////////////
// GETTERS

void LidDrivenCavity::getCoords(int* coords)
{
    coords[0] = this -> coords[0];
    coords[1] = this -> coords[1];
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

//////////////////////////////////////////////////////////////
// SOLVERS
void LidDrivenCavity::Initialise()
{
}

void LidDrivenCavity::Integrate()
{
}