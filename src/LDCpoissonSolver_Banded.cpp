#include <iostream>
#include <iomanip>
#include <cstring>
#include <mpi.h>

#include "LDCpoissonSolver_Banded.h"

#define F77NAME(x) x##_
extern "C" {
    double F77NAME(dpbtrf) (const char& UPLO, const int& n,
                            const int& kd, const double* AB,
                            const int& ldab, int& info);

    double F77NAME(dpbtrs) (const char& UPLO, const int& n,
                            const int& kd, const int& nRHS,
                            const double* AP, const int& ldab,
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

    LDCpoissonSolver_Banded::~LDCpoissonSolver_Banded(){
        delete[] A;
        delete[] b;
        }

//////////////////////////////////////////////////////////////
// SOLVERS

    void LDCpoissonSolver_Banded::Initialize(int& Nx, int& Ny, double* coeff){

        if(rank == 0){
            cout << "Initializing Poisson Solver" << endl << endl;
        }

        this -> Nx = Nx-2;
        this -> Ny = Ny-2;
        nNodes = (this -> Nx)*(this -> Ny); // Number of nodes in domain
        nCoeffs = nNodes * (this -> Ny+1);

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

    void LDCpoissonSolver_Banded::Build2DLaplace(){

        if(rank == 0){
            cout << "Building Coefficient matrix - Banded storage" <<  endl << endl;
        }

        // First diagonal entry
        int kd = Ny +1; // Number of superdiagonals

        A[Ny] = coeff[1];

        // Populate coefficient matrix
        for (int i=1; i<nNodes; i++) {

            // Diagonal entries
            A[(i+1)*kd - 1] = coeff[1];

            // Superdiagonal - y-direction
            if (i % Ny != 0){
                A[(i+1)*kd - 2] = coeff[0];
            }

            // Superdiagonal - x-direction
            if (i >= Ny - 1){
                A[(i+1)*kd] = coeff[2];
            }
        }
    }

    void LDCpoissonSolver_Banded::BuildRHS(double* v, double* s){

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

    void LDCpoissonSolver_Banded::SolvePoisson(double* v, double* s){

        // Update RHS vector
        this -> BuildRHS(v,s);

        // Solve linear system with LAPACK:
        // Packed storage, Symetric Positive definite matrix
        F77NAME(dpbtrs) ('U', nNodes, Ny, 1, A, Ny +1, b, nNodes, info);

        // Position solution back in streamfunction array
        int offset = Ny + 2;
        for(int i =0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &b[i*Ny], 1, &s[offset*(i+1) + 1], 1);
        }

    }

    void LDCpoissonSolver_Banded::Factor2DLaplace(){

        if (rank ==0 ){
            cout << "Computing coefficient matrix factor" << endl;
        }

        F77NAME(dpbtrf) ('U', nNodes, Ny, A, Ny+1, info);
    }

//////////////////////////////////////////////////////////////
// IO METHODS

    void LDCpoissonSolver_Banded::PrintCoeffMat(int rank){
        if (this -> rank == rank ){

            cout << "Nx=" << Nx << "  -  Ny=" << Ny  << "  -  nNodes=" << nNodes << "  -  nCoeffs=" << nCoeffs << endl << endl;

            for (int j = 0; j < Ny+1 ; j++){
                for (int i = 0; i < nNodes; i++){
                    cout << setw(8) << setprecision(3) << A[j + i*(Ny + 1) ];
                }
                cout << endl;
            }

        }
    }

    void LDCpoissonSolver_Banded::PrintRHS(int rank) {

        if (this -> rank == rank){

            for (int j=0; j < Ny + 1 ; j++){
                for (int i=0; i < Nx; i++){
                    cout << setw(8) << setprecision(3) << b[j + Ny*i] << "  ";
                }
                cout << endl;
            }

            cout << endl;
        }
    }
