#include <iostream>
#include <iomanip>
#include <cstring>
#include <mpi.h>

#include "LDCpoissonSolver_CGS.h"

#define F77NAME(x) x##_

// Definition of BLAS functions
extern "C" {
    // Matrix Vector Multiplication
    // Matrix stored in symmetric banded format
    double F77NAME(dsbmv) (const char& UPLO, const int& n,
                          const int& k, const double& alpha,
                          const double* A, const int& lda,
                          const double* x, const int& incx,
                          const double& beta, const double* y, 
                          const int& incy);

    // Y = X
    double F77NAME(dcopy) (const int& n,
                          const double* x, const int& incx,
                          const double* y, const int& incy);

    // Y = alpha*X + Y
    double F77NAME(daxpy) (const int& n, const double& alpha,
                          const double* x, const int& incx,
                          const double* y, const int& incy);

    // ddot = X_T * X
    double F77NAME(ddot) (const int& n, const double* x,
                         const int& incx, const double* y,
                         const int& incy);

    // Y = alpha * Y
    double F77NAME(dscal) (const int& n, const double& alpha,
                          const double* x, const int& incx);
}

//////////////////////////////////////////////////////////////
// CONSTRUCTORS

    LDCpoissonSolver_CGS::~LDCpoissonSolver_CGS(){
        delete[] A;
        delete[] Ap;
        delete[] p;
        delete[] r;
        delete[] u;
        delete[] matvecbuf_Nx;
        delete[] matvecbuf_Ny;
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
        // in BLAS banded storage format
        A  = new double [nCoeffs]{};
        Ap = new double [nNodes]{};
        p  = new double [nNodes]{};
        r  = new double [nNodes]{};
        u  = new double [nNodes]{};

        // Assign memory to parallel MatVec multiplication buffer arrays
        matvecbuf_Nx = new double [this -> Nx]{};
        matvecbuf_Ny = new double [this -> Ny]{};


        // Compute coefficients of the 2D FD Laplacian operator
        this -> coeff[0] = coeff[0];
        this -> coeff[1] = coeff[1];
        this -> coeff[2] = coeff[2];

        // Initialize 2D Laplace coefficient matrix
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

        int offset = Ny + 2;

        // Populate vorticity forcing values
        for (int i = 0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &v[offset*(i+1) + 1], 1, &r[i*Ny], 1);
        }
    }

    void LDCpoissonSolver_CGS::SolvePoisson(double* v, double* s){
        // Implements a parallel Conjugate Gradient solver

        // Declare local solver variables
        int offset = Ny + 2, k = 0;
        double dotR, dotR_prev, dotP;
        double alpha, beta;

        // Update forcing vector
        this -> BuildRHS(v,s);

        // Initialize u vector with internal streamfunction values
        for (int i = 0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &s[offset*(i+1) + 1], 1, &u[i*Ny], 1);
        }

        // Initialize r vector, adding contribution of neighbor nodes to mat-vec product
            F77NAME(dsbmv) ('U', nNodes, Ny, -1.0, A, Ny + 1, u, 1, 1.0, r, 1);

            // Neigbhor values x-direction
            F77NAME(daxpy) (Ny, -coeff[2], &s[1], 1, r, 1);
            F77NAME(daxpy) (Ny, -coeff[2], &s[offset*(Nx+1) + 1], 1, &r[Ny*(Nx-1)], 1);

            // Neighbor values y-direction
            F77NAME(daxpy) (Nx, -coeff[0], &s[Ny +2], (Ny + 2), r, Ny);
            F77NAME(daxpy) (Nx, -coeff[0], &s[2*(Ny +2) - 1], (Ny + 2), &r[Ny-1], Ny);

        F77NAME(dcopy) (nNodes, r, 1, p, 1);

        // Compute local/global r_prev-dot products
        dotR_prev = F77NAME(ddot) (nNodes, r, 1, r, 1);
        MPI_Allreduce(MPI_IN_PLACE, &dotR_prev, 1, MPI_DOUBLE, MPI_SUM, MPIcomm);

        /////////////// ITERATE
        while (true){

            // Update iteration count
            k++;

            // Reset and Update Ap vector
            fill_n(Ap,nNodes,0.0);
            InterfaceBroadcast(p);
            InterfaceGather(Ap);

            // Compute pTAp
            F77NAME(dsbmv) ('U', nNodes, Ny, 1.0, A, Ny + 1, p, 1, 1.0, Ap, 1);

            // Compute local/global p-dot products
            dotP = F77NAME(ddot) (nNodes, p, 1, Ap, 1);
            MPI_Allreduce(MPI_IN_PLACE, &dotP, 1, MPI_DOUBLE, MPI_SUM, MPIcomm);

            // Update alpha
            alpha = dotR_prev / dotP;

            // Update u
            F77NAME(daxpy) (nNodes, alpha, p, 1, u, 1);

            // Update r
            F77NAME(daxpy) (nNodes, -alpha, Ap, 1, r, 1);

            // Compute local/global r-dot products
            dotR = F77NAME(ddot) (nNodes, r, 1, r, 1);
            MPI_Allreduce(MPI_IN_PLACE, &dotR, 1, MPI_DOUBLE, MPI_SUM, MPIcomm);

            // Check convergence and max iteration criteria
            int iter_max = 1000; // Adjust as needed
            double res_crit = 0.0000001; // Adjust as needed

            if (k > iter_max || dotR < res_crit){
                break;
            }

            // Update convergence criteria
            beta = dotR / dotR_prev;
            dotR_prev = dotR;

            // Update p vector
            F77NAME(dscal) (nNodes, beta, p, 1);
            F77NAME(daxpy) (nNodes, 1, r, 1, p, 1);

        } 

        // Position solution back in streamfunction array
        for(int i =0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &u[i*Ny], 1, &s[offset*(i+1) + 1], 1);
        }

        if (rank == 0){
            cout << "Solver converged to tolerance in " << k << " iterations" << endl;
        }
    }

//////////////////////////////////////////////////////////////
// MPI INTERFACE MANAGEMENT

    void LDCpoissonSolver_CGS::InterfaceBroadcast(double* field){

        //Sequentially send interface values to neighbors in all directions

        // Neighbor below
        if (rankShift[0] != -2){
            LDCpoissonSolver_CGS::InterfaceSend(Nx, field, matvecbuf_Nx, Ny, rankShift[0], rank, MPIcomm);
        }

        // Neighbor above
        if (rankShift[1] != -2){
            LDCpoissonSolver_CGS::InterfaceSend(Nx, &field[Ny-1], matvecbuf_Nx, Ny, rankShift[1], rank, MPIcomm);
        }

        // Neighbor leftward
        if (rankShift[2] != -2){
            LDCpoissonSolver_CGS::InterfaceSend(Ny, field, matvecbuf_Ny, 1, rankShift[2], rank, MPIcomm);
        }

        // Neighbor rightward
        if (rankShift[3] != -2){
            LDCpoissonSolver_CGS::InterfaceSend(Ny, &field[Ny*(Nx-1)], matvecbuf_Ny, 1, rankShift[3], rank, MPIcomm);
        }
    }

    void LDCpoissonSolver_CGS::InterfaceGather(double* field){

        //Sequentially send interface values to neighbors in all directions

        //Neighbor above
        if (rankShift[1] != -2){
            LDCpoissonSolver_CGS::InterfaceRecv(Nx, coeff[0], &field[Ny-1], matvecbuf_Nx, Ny, rankShift[1], rankShift[1], MPIcomm);
        }

        //Neighbor below
        if (rankShift[0] != -2){
            LDCpoissonSolver_CGS::InterfaceRecv(Nx, coeff[0], field, matvecbuf_Nx, Ny, rankShift[0], rankShift[0], MPIcomm);
        }

        //Neighbor rightward
        if (rankShift[3] != -2){
            LDCpoissonSolver_CGS::InterfaceRecv(Ny, coeff[2], &field[Ny*(Nx-1)], matvecbuf_Ny, 1, rankShift[3], rankShift[3], MPIcomm);
        }

        //Neighbor leftward
        if (rankShift[2] != -2){
            LDCpoissonSolver_CGS::InterfaceRecv(Ny, coeff[2], field, matvecbuf_Ny, 1, rankShift[2], rankShift[2], MPIcomm);
        }
    }

    // Copy into buffer and send to target
    void LDCpoissonSolver_CGS::InterfaceSend(int& count, double* field, double* buff, int disp, int& dest, int& tag, MPI_Comm MPIcomm){
        F77NAME(dcopy) (count, field, disp, buff, 1);
        MPI_Send(buff, count, MPI_DOUBLE, dest, tag, MPIcomm);
    }

    // Receive into buffer then scale and add to target field
    void LDCpoissonSolver_CGS::InterfaceRecv(int& count, double& alpha, double* field, double* buff, int disp, int& source, int& tag, MPI_Comm MPIcomm){
        MPI_Recv(buff, count, MPI_DOUBLE, source, tag, MPIcomm, MPI_STATUS_IGNORE);
        F77NAME(daxpy) (count,alpha, buff, 1, field, disp);
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
