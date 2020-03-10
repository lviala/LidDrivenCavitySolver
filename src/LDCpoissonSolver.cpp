#include <iostream>
#include <iomanip>
#include <cstring>
#include <mpi.h>

#include "LDCpoissonSolver.h"

#define F77NAME(x) x##_
extern "C" {
    double F77NAME(dcopy) (const int& n,
                          const double *x, const int& incx,
                          const double *y, const int& incy);

    double F77NAME(daxpy) (const int& n, const double& alpha,
                          const double *x, const int& incx,
                          const double *y, const int& incy);
}

//////////////////////////////////////////////////////////////
// CONSTRUCTORS
    LDCpoissonSolver::LDCpoissonSolver(int rank){
        this -> rank = rank;
    }

//////////////////////////////////////////////////////////////
// SOLVERS

    void LDCpoissonSolver::BuildRHS(double* v, double* s){

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


//////////////////////////////////////////////////////////////
// IO METHODS

    void LDCpoissonSolver::PrintRHS(int rank) {

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