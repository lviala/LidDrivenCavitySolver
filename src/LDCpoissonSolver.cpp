#include <iostream>
#include <iomanip>
#include <cstring>
#include <mpi.h>

#include "LDCpoissonSolver.h"

//////////////////////////////////////////////////////////////
// SETTERS
LDCpoissonSolver::LDCpoissonSolver(int rank){
    this -> rank = rank;
}

LDCpoissonSolver::~LDCpoissonSolver(){
    delete A;
}

//////////////////////////////////////////////////////////////
// SOLVERS

void LDCpoissonSolver::Initialize(int& Nx, int& Ny, double& dx, double& dy){

    this -> Nx = Nx-2;
    this -> Ny = Ny-2;
    nNodes = (this -> Nx)*(this -> Ny); // Number of nodes in domain
    nCoeffs = nNodes*(nNodes+1)/2;
    
    // Assign memory to coefficient Matrix A
    // in LAPACK packed storage format
    A = new double [nCoeffs];
    fill_n(A,nCoeffs,0.0);

    // Compute coefficients of the 2D FD Laplacian operator
    coeff[0] = -1.0/(dy*dy);
    coeff[2] = -1.0/(dx*dx);
    coeff[1] = -(coeff[0] + coeff[2]);

    this -> Build2DLaplace();
    
}

void LDCpoissonSolver::Build2DLaplace(){

    // Index of diagonal entry
    int idDiag;

    for (int i=0; i<nNodes; i++) {

        idDiag = (i+1)*(i+2)/2 -1 ;
        //  Populate diagonal entries
        A[idDiag] = coeff[1];

        // Populate y-direction neighbor
        if (i % Ny != 0){
            A[idDiag - 1] = coeff[0];
        }

        if (i >= Ny) {
            // Populate x-direction neighbor
            A[idDiag - Ny] = coeff[2];
        }

    }
    if (rank == 0){

        cout << "Nx=" << Nx << "  -  Ny=" << Ny  << "  -  nNodes=" << nNodes << "  -  nCoeffs=" << nCoeffs << endl << endl;

        
        for (int i = 0; i < nNodes; i++){
            for (int j = 0; j < nNodes; j++){
                
                if (i <= j){
                    cout << setw(6) << A[(i) + j*(j+1)/2];
                }
                else {
                    cout << setw(6) << "";
                }
            }
            cout << endl;
        }
        
    }
    
}