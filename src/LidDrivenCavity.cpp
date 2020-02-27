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

void LidDrivenCavity::getStreamFunction(double* s)
{
    s = this -> s;
}
void LidDrivenCavity::getVorticity(double* v)
{
    v = this -> v;
}
void LidDrivenCavity::getReynoldsNumber(double& Re)
{
    Re = this -> Re;
}

//////////////////////////////////////////////////////////////
// SOLVERS
void LidDrivenCavity::Initialise()
{
    // Compute the necessary overlap required between subdomains
    // Based on neighbor of cartesian subgrid
    // rankShift[i] = -2 indicates MPI_PROC_NULL

    if (rankShift[0] != -2){
        Ny += 1;
    }
    if (rankShift[1] != -2){
        Ny += 1;
    }
    if (rankShift[2] != -2){
        Nx += 1;
    }
    if (rankShift[3] != -2){
        Nx += 1;
        }

    // Allocate memory to subgrid vorticity and streamfunction fields
    this -> v = new double [Nx*Ny];
    memset(this -> v , 0 , Ny * Nx );
    this -> s = new double [Nx*Ny];
    memset(this -> s , 0 , Ny * Nx);

}

void LidDrivenCavity::Integrate()
{
}

//////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

void LidDrivenCavity::PrintArray(const char* varStr, int rank) {

    if (this -> rank == rank){
        double* toPrint = nullptr;
        
        if (strcmp(varStr,"s") == 0){
            toPrint = this -> s;
        }
        else if (strcmp(varStr,"v") == 0){
            toPrint = this -> v;
        }

        for (int j=0; j < Ny; j++){
            for (int i=0; i < Nx; i++){
                cout << toPrint[j + (Nx-1)*i] << "  ";
            }
            cout << endl;
        }
        
        cout << endl;
    }
}

void LidDrivenCavity::LDCStatus(int rank){
    if (this -> rank == rank){
        cout << "My rank is: " << this -> rank << endl; 
        cout << "My Coordinates are: (" << this -> coords[0] << " , " << this -> coords[1] << ")" << endl;
        cout << "My Neighbors are: Down=" << rankShift[0] << "  -- Up =" << rankShift[1] << "  -- Left =" << rankShift[2] << "  -- Right =" << rankShift[3] << endl;
        cout << "Nx=" << this -> Nx << " -- Ny=" << this -> Ny << endl << endl;
    }
}