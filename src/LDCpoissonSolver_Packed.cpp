#include <iostream>
#include <iomanip>
#include <cstring>
#include <mpi.h>

#include "LDCpoissonSolver_Packed.h"

#define F77NAME(x) x##_
extern "C" {
    double F77NAME(dpptrf) (const char& UPLO, const int& n,
                            const double* AP, int& info);

    double F77NAME(dpptrs) (const char& UPLO, const int& n,
                            const int& nRHS, const double* AP, 
                            const double* b, const int& ldb,
                            int& info);

    double F77NAME(dcopy) (const int& n,
                          const double *x, const int& incx,
                          const double *y, const int& incy);

    double F77NAME(daxpy) (const int& n, const double& alpha,
                          const double *x, const int& incx,
                          const double *y, const int& incy);
}

//////////////////////////////////////////////////////////////
// CONSTRUCTORS
    // LDCpoissonSolver_Packed::LDCpoissonSolver_Packed(int rank){
    // }

    LDCpoissonSolver_Packed::~LDCpoissonSolver_Packed(){
        delete[] A;
        delete[] b;
        }

//////////////////////////////////////////////////////////////
// SOLVERS

    void LDCpoissonSolver_Packed::Initialize(int& Nx, int& Ny, double* coeff){

        if(rank == 0){
            cout << "Initializing Poisson Solver" << endl << endl;
        }

        this -> Nx = Nx-2;
        this -> Ny = Ny-2;
        nNodes = (this -> Nx)*(this -> Ny); // Number of nodes in domain
        nCoeffs = nNodes*(nNodes+1)/2;
        
        // Assign memory to coefficient Matrix A
        // in LAPACK packed storage format
        A = new double [nCoeffs]{};
        b = new double [nNodes]{};

        // Compute coefficients of the 2D FD Laplacian operator
        this -> coeff[0] = coeff[0];
        this -> coeff[1] = coeff[1];
        this -> coeff[2] = coeff[2];

        // Initialize 2D Laplace coefficient matrix and factor
        this -> Build2DLaplace();
        this -> Factor2DLaplace();
        
    }

    void LDCpoissonSolver_Packed::Build2DLaplace(){

        if(rank == 0){
            cout << "Building Coefficient matrix" <<  endl << endl;
        }

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

    void LDCpoissonSolver_Packed::BuildRHS(double* v, double* s){

        // Reset RHS vector
        fill_n(b,nNodes,0.0);
        int offset = Ny + 2;

        // Populate vorticity forcing values
        for (int i = 0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &v[offset*(i+1) + 1], 1, &b[i*Ny], 1);
        }

        // Populate streamfunction BC values x-direction
        F77NAME(daxpy) (Ny, -coeff[2], &s[1], 1, b, 1);
        F77NAME(daxpy) (Ny, -coeff[2], &s[offset*(Nx+1) + 1], 1, &b[Ny*(Nx-1)], 1);

        // Populate streamfunction BC values y-direction
        F77NAME(daxpy) (Nx, -coeff[0], &s[Ny +2], (Ny + 2), b, Ny);
        F77NAME(daxpy) (Nx, -coeff[0], &s[2*(Ny +2) - 1], (Ny + 2), &b[Ny-1], Ny);

    }

    void LDCpoissonSolver_Packed::SolvePoisson(double* v, double* s){

        // Update RHS vector
        this -> BuildRHS(v,s);

        // Solve linear system with LAPACK:
        // Packed storage, Symetric Positive definite matrix
        F77NAME(dpptrs) ('U', nNodes, 1, A, b, nNodes, info);

        // Position solution back in streamfunction array
        int offset = Ny + 2;
        for(int i =0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &b[i*Ny], 1, &s[offset*(i+1) + 1], 1);
        }

    }

    void LDCpoissonSolver_Packed::Factor2DLaplace(){
        
        if (rank ==0 ){
            cout << "Computing coefficient matrix factor" << endl;
        }

        F77NAME(dpptrf) ('U', nNodes, A, info);
        
    }

//////////////////////////////////////////////////////////////
// IO METHODS

    void LDCpoissonSolver_Packed::PrintCoeffMat(int rank){
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

    void LDCpoissonSolver_Packed::PrintRHS(int rank) {

        if (this -> rank == rank){
            
            for (int j=0; j < Ny; j++){
                for (int i=0; i < Nx; i++){
                    cout << setw(8) << setprecision(3) << b[j + Ny*i] << "  ";
                }
                cout << endl;
            }
            
            cout << endl;
        }
    }