#include <iostream>
#include <iomanip>
#include <cstring>
#include <mpi.h>

#include "LDCpoissonSolver_CGS.h"

#define F77NAME(x) x##_
extern "C" {
    double F77NAME(dsbmv) (const char& UPLO, const int& n,
                          const int& k, const double alpha,
                          const double* A, const int& lda,
                          const double* x, const int& incx,
                          const double& beta, const double* y, 
                          const int& incy);

    double F77NAME(dcopy) (const int& n,
                          const double *x, const int& incx,
                          const double *y, const int& incy);

    double F77NAME(daxpy) (const int& n, const double& alpha,
                          const double *x, const int& incx,
                          const double *y, const int& incy);

    double F77NAME(ddot) (const int& n, const double* x,
                         const int& incx, const double* y,
                         const int& incy);
}

//////////////////////////////////////////////////////////////
// CONSTRUCTORS

    LDCpoissonSolver_CGS::~LDCpoissonSolver_CGS(){
        delete[] A;
        delete[] b;
        }

//////////////////////////////////////////////////////////////
// SOLVERS

    void LDCpoissonSolver_CGS::Initialize(int& Nx, int& Ny, double* coeff){

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
        p = new double [nNodes]{};
        r = new double [nNodes]{};
        u = new double [nNodes]{};

        // Compute coefficients of the 2D FD Laplacian operator
        this -> coeff[0] = coeff[0];
        this -> coeff[1] = coeff[1];
        this -> coeff[2] = coeff[2];

        // Initialize 2D Laplace coefficient matrix and factor
        this -> Build2DLaplace();
    }

    void LDCpoissonSolver_CGS::Build2DLaplace(){

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

    void LDCpoissonSolver_CGS::BuildRHS(double* v, double* s){

        // Reset RHS vector
        fill_n(b,nNodes,0.0);
        int offset = Ny + 2;

        // Populate vorticity forcing values
        for (int i = 0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &v[offset*(i+1) + 1], 1, &r[i*Ny], 1);
        }
    }

    void LDCpoissonSolver_CGS::SolvePoisson(double* v, double* s){

        int offset = Ny + 2, k = 0;

        double dotR, dotR_prev, dotP;
        double alpha, beta;

        // Update forcing vector
        this -> BuildRHS(v,s);

        // Initialize u vector with internal streamfunction values
        for (int i = 0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &s[offset*(i+1) + 1], 1, &u[i*Ny], 1);
        }

        // Initialize r vector, p vector
        F77NAME(dsbmv) ('U', nNodes, Ny, -1.0, A, nNodes, u, 1, -1.0, r, 1);
        F77NAME(dcopy) (nNodes, r, 1, p, 1);

        // Compute local/global r_prev-dot products
        dotR_prev = F77NAME(ddot) (nNodes, r, 1, r, 1);
        MPI_Allreduce(MPI_IN_PLACE, &dotR_prev, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


        /////////////// ITERATE

        while (true){
            // Update p vector

            // Compute local/global p-dot products
            dotP = F77NAME(ddot) (nNodes, p, 1, p, 1);
            MPI_Allreduce(MPI_IN_PLACE, &dotP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            // Update alpha
            alpha = dotR_prev / dotP;

            // Update u
            F77NAME(daxpy) (nNodes, alpha, p, 1, u, 1);

            // Update r

            // Compute local/global r-dot products
            dotR = F77NAME(ddot) (nNodes, r, 1, r, 1);
            MPI_Allreduce(MPI_IN_PLACE, &dotR, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            if (k < 1000 || dotR < 1e-3){
                break;
            }

            // Update convergence criteria
            beta = dotR / dotR_prev;
            dotR_prev = dotR;

            // Update p vector
            F77NAME(daxpy) (nNodes, )

            k++;

        } 
        



        
        // Populate streamfunction BC values x-direction
        F77NAME(daxpy) (Ny, -coeff[2], &s[1], 1, b, 1);
        F77NAME(daxpy) (Ny, -coeff[2], &s[offset*(Nx+1) + 1], 1, &b[Ny*(Nx-1)], 1);

        // Populate streamfunction BC values y-direction
        F77NAME(daxpy) (Nx, -coeff[0], &s[Ny +2], (Ny + 2), b, Ny);
        F77NAME(daxpy) (Nx, -coeff[0], &s[2*(Ny +2) - 1], (Ny + 2), &b[Ny-1], Ny);

        // Position solution back in streamfunction array
        int offset = Ny + 2;
        for(int i =0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &b[i*Ny], 1, &s[offset*(i+1) + 1], 1);
        }
    }


//////////////////////////////////////////////////////////////
// IO METHODS

    void LDCpoissonSolver_CGS::PrintCoeffMat(int rank){
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

    void LDCpoissonSolver_CGS::PrintRHS(int rank) {

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
