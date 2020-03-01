#include "LidDrivenCavity.h"
#include "LDCpoissonSolver.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <mpi.h>

#define F77NAME(x) x##_
extern "C" {
    double F77NAME(dcopy) (const int& n,
                          const double *x, const int& incx,
                          const double *y, const int& incy);
}

//////////////////////////////////////////////////////////////
// CONSTRUCTORS

LidDrivenCavity::LidDrivenCavity(MPI_Comm MPIcomm, int rank, int* rankShift, int* coords, int* gridSize, double dt, double dx, double dy, double T, double Re)
{   
    this -> MPIcomm = MPIcomm;
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
    this -> dx = dx;
    this -> dy = dy;
    this -> T = T;
    this -> Re = Re;
}

LidDrivenCavity::~LidDrivenCavity()
{
    delete poissonSolver;
    delete v;
    delete s;
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
    // Based on neighbors of cartesian subgrid
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

    // Allocate memory to fields and buffers
    this -> v = new double [Nx*Ny];
    fill_n(this -> v , Ny * Nx, 0 );

    this -> s = new double [Nx*Ny];
    fill_n(this -> v , Ny * Nx, 0 );


    this -> bufNx = new double [Nx];
    fill_n(this -> bufNx , Nx, 0);
    this -> bufNy = new double [Ny];
    fill_n(this -> bufNy , Ny, 0);

    poissonSolver = new LDCpoissonSolver(rank);
    poissonSolver -> Initialize(Nx, Ny, dx, dy);

}

void LidDrivenCavity::UpdateGlobalBcs(){
    // Determine if subgrid is on the edge of the global domain
    // and set appropriate BCs on vorticity field.
    // Recall direction convention of array: 
    //      Column major with [+ve y-direction -> down , +ve x_direction -> right]

    // Domain on bottom of cavity if rankshift[0] = -2
    if (rankShift[0] == -2){
        for (int i = 0; i < this -> Nx ; i++){
            this -> v[i*Ny] = (2/(dy*dy))*(s[i*(Ny)] - s[1 + i*(Ny)]);
        }
    }

    // Domain on top of cavity if rankshift[1] = -2
    if (rankShift[1] == -2){
        for (int i = 0; i < this -> Nx ; i++){
            this -> v[(Ny-1) + i*Ny] = (2/(dy*dy))*(s[(Ny-1) + i*(Ny)] - s[(Ny-2) + i*(Ny)]) - 2.0*U/dy;
        }
    }

    // Domain on left of cavity if rankshift[2] = -2
    if (rankShift[2] == -2){
        for (int j = 0; j < this -> Ny; j++){
            this -> v[j] = (2/(dx*dx))*(s[j] - s[j + Ny]);
        }
    }

    // Domain on right of cavity if rankshift[3] = -2
    if (rankShift[3] == -2){
        for (int j = 0; j < this -> Ny; j++){
            this -> v[j + (Nx - 1)*Ny] = (2/(dx*dx))*(s[j+ (Nx - 1)*Ny] - s[j + (Nx - 2)*Ny]);
        }
    }
}

void LidDrivenCavity::Integrate()
{
    InterfaceBroadcast(s);
    InterfaceGather(s);
}

//////////////////////////////////////////////////////////////
// MPI INTERFACE MANAGEMENT

void LidDrivenCavity::InterfaceBroadcast(double* field){
    
    //Sequentially send interface values to neighbors in all directions
    
    // Neighbor below
    if (rankShift[0] != -2){
        LidDrivenCavity::InterfaceSend(Nx, field, bufNx, Ny, rankShift[0],rank,MPIcomm);
    }

    // Neighbor above
    if (rankShift[1] != -2){
        LidDrivenCavity::InterfaceSend(Nx, &field[Ny-1], bufNx, Ny, rankShift[1],rank,MPIcomm);
    }

    // Neighbor leftward
    if (rankShift[2] != -2){
        LidDrivenCavity::InterfaceSend(Ny, field, bufNy, 1, rankShift[2],rank,MPIcomm);
    }

    // Neighbor rightward
    if (rankShift[3] != -2){
        LidDrivenCavity::InterfaceSend(Ny, &field[Ny*(Nx-1)], bufNy, 1, rankShift[3],rank,MPIcomm);
    }

}

void LidDrivenCavity::InterfaceGather(double* field){
    
    //Sequentially send interface values to neighbors in all directions

    //Neighbor above
    if (rankShift[1] != -2){
        LidDrivenCavity::InterfaceRecv(Nx, &field[Ny-1], bufNx, Ny, rankShift[1], rankShift[1], MPIcomm);
    }

    //Neighbor below
    if (rankShift[0] != -2){
        LidDrivenCavity::InterfaceRecv(Nx, field, bufNx, Ny, rankShift[0], rankShift[0], MPIcomm);
    }

    //Neighbor rightward
    if (rankShift[3] != -2){
        LidDrivenCavity::InterfaceRecv(Ny, &field[Ny*(Nx-1)], bufNy, 1, rankShift[3], rankShift[3], MPIcomm);
    }

    //Neighbor leftward
    if (rankShift[2] != -2){
        LidDrivenCavity::InterfaceRecv(Ny, field, bufNy, 1, rankShift[2], rankShift[2], MPIcomm);
    }
}


void LidDrivenCavity::InterfaceSend(int& count, double* field, double* buff, int disp, int& dest, int& tag, MPI_Comm MPIcomm){
    F77NAME(dcopy) (count, field, disp, buff, 1);
    MPI_Send(buff, count, MPI_DOUBLE, dest, tag, MPIcomm);
}

void LidDrivenCavity::InterfaceRecv(int& count, double* field, double* buff, int disp, int& source, int& tag, MPI_Comm MPIcomm){
    MPI_Recv(buff, count, MPI_DOUBLE, source, tag, MPIcomm, MPI_STATUS_IGNORE);
    F77NAME(dcopy) (count, buff, 1, field, disp);
}

//////////////////////////////////////////////////////////////
// IO FUNCTIONS

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
                cout << setw(4) << toPrint[j + Ny*i] << "  ";
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