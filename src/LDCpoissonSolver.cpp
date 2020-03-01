#include <iostream>
#include <iomanip>
#include <cstring>
#include <mpi.h>

#include "LDCpoissonSolver.h"

#define F77NAME(x) x##_
extern "C" {
    double F77NAME(dpptrf) (const char& UPLO, const int& n,
                            const double* AP, int& info);

    double F77NAME(dpptrs) (const char UPLO, const int& n,
                            const int& nRHS, const double* AP, 
                            const double* b, const int& ldb,
                            int& info);
}

//////////////////////////////////////////////////////////////
// CONSTRUCTORS
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
    A = new double [nCoeffs]{};

    // Compute coefficients of the 2D FD Laplacian operator
    coeff[0] = -1.0/(dy*dy);
    coeff[2] = -1.0/(dx*dx);
    coeff[1] = -2.0*(coeff[0] + coeff[2]);

    this -> Build2DLaplace();
    this -> PrintTRIUArray(0);
    this -> Factor2DLaplace();
    this -> PrintTRIUArray(0);
    
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
}

void LDCpoissonSolver::Factor2DLaplace(){
    if (rank == 0){
        cout << nNodes << endl;
        F77NAME(dpptrf) ('U', nNodes, A, info);
    }
}

void LDCpoissonSolver::PrintTRIUArray(int rank){
    if (this -> rank == rank ){

        cout << "Nx=" << Nx << "  -  Ny=" << Ny  << "  -  nNodes=" << nNodes << "  -  nCoeffs=" << nCoeffs << endl << endl;

        
        for (int i = 0; i < nNodes; i++){
            for (int j = 0; j < nNodes; j++){
                
                if (i <= j){
                    cout << setw(8) << setprecision(3) <<  A[(i) + j*(j+1)/2];
                }
                else {
                    cout << setw(8) << "";
                }
            }
            cout << endl;
        }
    
    }
}