#include "LidDrivenCavity.h"
#include "LDCpoissonSolver_Banded.h"
#include "LDCpoissonSolver_CGS.h"
#include <iostream>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <mpi.h>

#define F77NAME(x) x##_
// Definition of BLAS functions
extern "C" {
    // Y = X
    double F77NAME(dcopy) (const int& n,
                          const double *x, const int& incx,
                          const double *y, const int& incy);

    // Norm2 of vector
    double F77NAME(dnrm2) (const int& n, const double* x, const int& incx);

    // Y = alpha*X + Y
    double F77NAME(daxpy) (const int& n, const double& alpha,
                          const double *x, const int& incx,
                          const double *y, const int& incy);
}

//////////////////////////////////////////////////////////////
// CONSTRUCTORS

    LidDrivenCavity::LidDrivenCavity(MPI_Comm MPIcomm, int rank, int* rankShift, int* coords, int* gridSize, double dt, double dx, double dy, double* subPos, double T, double Re)
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

        this -> dx = dx;
        this -> dy = dy;

        this -> subPos[0] = subPos[0];
        this -> subPos[1] = subPos[1];

        this -> T = T;
        this -> Re = Re;
        this -> dt = dt;

    }

    LidDrivenCavity::~LidDrivenCavity()
    {
        // delete[] v;
        // delete[] v_new;
        // delete[] s;
        delete[] velU;
        delete[] velV;
        delete[] bufNx;
        delete[] bufNy;
        delete poissonSolver;
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

        // Compute coefficients of the 2D FD Laplacian operator
        coeff[0] = -1.0/(dy*dy);
        coeff[2] = -1.0/(dx*dx);
        coeff[1] = -2.0*(coeff[0] + coeff[2]);

        // Allocate memory to fields and buffers
        this -> s       = new double [Nx*Ny]{};
        this -> v       = new double [Nx*Ny]{};
        this -> v_new   = new double [Nx*Ny]{};
        this -> velU    = new double [Nx*Ny]{};
        this -> velV    = new double [Nx*Ny]{};
        this -> bufNx   = new double [Nx]{};
        this -> bufNy   = new double [Ny]{};

        // Set lid top velocity
        if (rankShift[1] == -2){
            for (int i = 1; i < Nx ; i++){
                velU[i*Ny - 1] = this -> U;
            }
        }

        // Initial update of global BCs
        this -> UpdateGlobalBcs();

        // Initialize poisson solver object
        // Builds coefficient matrix and other necessary variables

        poissonSolver = new LDCpoissonSolver_CGS(rank, rankShift, MPIcomm);
        poissonSolver -> Initialize(Nx, Ny, coeff);

    }

    void LidDrivenCavity::UpdateGlobalBcs(){
        // Determine if subgrid is on the edge of the global domain
        // and set appropriate BCs on vorticity field.
        // Recall direction convention of array:
        //      Column major with [+ve y-direction -> down , +ve x_direction -> right]

        // Domain on bottom of cavity if rankshift[0] = -2
        if (rankShift[0] == -2){
            for (int i = 0; i < this -> Nx ; i++){
                v[i*Ny] = (2.0/(dy*dy))*(s[i*(Ny)] - s[1 + i*(Ny)]);
            }
        }

        // Domain on top of cavity if rankshift[1] = -2
        if (rankShift[1] == -2){
            for (int i = 0; i < this -> Nx ; i++){
                v[(i+1)*Ny - 1] = (2.0/(dy*dy))*(s[(i+1)*Ny - 1] - s[(i+1)*Ny - 2]) - 2.0*U/dy;
            }
        }

        // Domain on left of cavity if rankshift[2] = -2
        if (rankShift[2] == -2){
            for (int j = 0; j < this -> Ny; j++){
                v[j] = (2.0/(dx*dx))*(s[j] - s[j + Ny]);
            }
        }

        // Domain on right of cavity if rankshift[3] = -2
        if (rankShift[3] == -2){
            for (int j = 0; j < this -> Ny; j++){
                v[j + (Nx - 1)*Ny] = (2.0/(dx*dx))*(s[j+ (Nx - 1)*Ny] - s[j + (Nx - 2)*Ny]);
            }
        }
    }

    void LidDrivenCavity::Integrate()
    {
        // Copy the values of v to v_new to preserve BCs
        F77NAME(dcopy)(Ny*Nx, v, 1, v_new, 1);

        // Overwrite the values of v_new with the laplacian of v
        FDLalplacianOperator(dt/Re, v, v_new);
        // Add the contribution of advection to v_new
        FDAdvectionOperator(dt, s, v, v_new);

        // Add the value of v to v_new
        F77NAME(daxpy) (Nx*Ny -2, 1.0, &v[1], 1, &v_new[1], 1);

        // Copy the values of v_new to v to preserve BCs
        F77NAME(dcopy)(Ny*Nx, v_new, 1, v, 1);
    }

    void LidDrivenCavity::Solve(){
        double t = 0.0;

        do{
            // Print solver status to console
            if(rank == 0){
                cout << "t = " << t << " of T = " << T << endl;
            }

            // Update vorticity field
            FDLalplacianOperator(-1.0, this -> s, this -> v);

            // Update interface values of the vorticity field
            InterfaceBroadcast(v);
            InterfaceGather(v);

            // Forward time step
            Integrate();

            // Update interface values of the vorticity field
            InterfaceBroadcast(v);
            InterfaceGather(v);

            // Solve the poisson problem to update the streamfunction field
            poissonSolver -> SolvePoisson(this -> v, this -> s);

            // Update interface values of the streamfunction field
            InterfaceBroadcast(s);
            InterfaceGather(s);

            // Update global boundary conditions
            UpdateGlobalBcs();

            // Update time
            t += dt;
        }
        while (t < T);

        // On solver completion, compute U and V velocity fields
        FDCurlOperator(-1, 1, s, velV, velU);
    }

//////////////////////////////////////////////////////////////
// FINITE DIFFERENCE OPERATORS

    void LidDrivenCavity::FDLalplacianOperator(const double& alpha, double* x, double* y){
        // Overwrites the values of array y with the laplacian of array x
        // Multiplied by scalar alpha

        for (int i=1; i < Nx-1; i++){
            for (int j=1; j< Ny-1; j++){
                y[ j + Ny*i ] = alpha*(( x[(j+1) + Ny*i] - 2.0*x[j + Ny*i] + x[(j-1) + Ny*i] )/(dy*dy) +
                                ( x[j + Ny*(i+1)] - 2.0*x[j + Ny*i] + x[j + Ny*(i-1)] )/(dx*dx));
            }
        }
    }

    void LidDrivenCavity::FDAdvectionOperator(const double& alpha, double* s, double* v, double* v_new){
        // Adds the advection component of the momentum equation to the values of array y

        for (int i=1; i < Nx-1; i++ ){
            for (int j=1; j< Ny-1; j++){
                v_new[j + Ny*i] += alpha * ((0.5/dy)*(v[(j+1) + Ny*i] - v[(j-1) + Ny*i]) * (0.5/dx)*(s[j + Ny*(i+1)] - s[j + Ny*(i-1)]) -
                                    (0.5/dy)*(s[(j+1) + Ny*i] - s[(j-1) + Ny*i]) * (0.5/dx)*(v[j + Ny*(i+1)] - v[j + Ny*(i-1)]));
            }
        }
    }

     void LidDrivenCavity::FDCurlOperator(double alpha_x, double alpha_y, double* f, double* df_dx, double* df_dy){
        for (int i=1; i < Nx-1; i++ ){
            for (int j=1; j< Ny-1; j++){
                df_dx[j + Ny*i] = alpha_x*(0.5/dx)*(f[j + Ny*(i+1)] - f[j + Ny*(i-1)]);
                df_dy[j + Ny*i] = alpha_y*(0.5/dy)*(f[(j+1) + Ny*i] - f[(j-1) + Ny*i]);
            }
        }
     }

//////////////////////////////////////////////////////////////
// MPI INTERFACE MANAGEMENT

    void LidDrivenCavity::InterfaceBroadcast(double* field){

        //Sequentially send interface values to neighbors in all directions

        // Neighbor below
        if (rankShift[0] != -2){
            LidDrivenCavity::InterfaceSend(Nx, &field[1], bufNx, Ny, rankShift[0],rank,MPIcomm);
        }

        // Neighbor above
        if (rankShift[1] != -2){
            LidDrivenCavity::InterfaceSend(Nx, &field[Ny-2], bufNx, Ny, rankShift[1],rank,MPIcomm);
        }

        // Neighbor leftward
        if (rankShift[2] != -2){
            LidDrivenCavity::InterfaceSend(Ny, &field[Ny], bufNy, 1, rankShift[2],rank,MPIcomm);
        }

        // Neighbor rightward
        if (rankShift[3] != -2){
            LidDrivenCavity::InterfaceSend(Ny, &field[Ny*(Nx-2)], bufNy, 1, rankShift[3],rank,MPIcomm);
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

    // Copy into buffer and send to target
    void LidDrivenCavity::InterfaceSend(int& count, double* field, double* buff, int disp, int& dest, int& tag, MPI_Comm MPIcomm){
        F77NAME(dcopy) (count, field, disp, buff, 1);
        MPI_Send(buff, count, MPI_DOUBLE, dest, tag, MPIcomm);
    }

    // Receive into buffer and copy into target field
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
                cout << "My Rank: " << this -> rank << " -- Streamfunction" << endl << endl;
                toPrint = this -> s;
            }
            else if (strcmp(varStr,"v") == 0){
                cout << "My Rank: " << this -> rank << " -- Vorticity" << endl << endl;
                toPrint = this -> v;
            }

            for (int j=0; j < Ny; j++){
                for (int i=0; i < Nx; i++){
                    cout << setw(8) << setprecision(3) << toPrint[j + Ny*i] << "  ";
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
            cout << "Nx=" << this -> Nx << " -- Ny=" << this -> Ny << " -- dx=" << this -> dx << " -- dy=" << this -> dy
                        << " -- Posx=" << this -> subPos[1] << " -- Posy=" << this -> subPos[0] << endl << endl;
        }
    }

    void LidDrivenCavity::LDCPrintSolution2File(string filename, const double& wallTime){
        // Sequentially write solution to csv file
        // Every process opens, writes into, and closes file
        // Only one process can have the file open at one time

        // Number of processes
        int size;
        MPI_Comm_size(MPIcomm, &size);

        // Index shifter to take into account interface nodes
        int yShift_Start = 0, yShift_End = 0;
        int xShift_Start = 0, xShift_End = 0;

        if (rankShift[0] != -2) yShift_Start++;
        if (rankShift[1] != -2) yShift_End++;
        if (rankShift[2] != -2) xShift_Start++;
        if (rankShift[3] != -2) xShift_End++;

        for (int k = 0; k < size; k++){
            if (k == 0 && rank == 0){ // Only for rank 0 and first write operation
                // Open file and set as output only and overwrite
                ofstream vOut(filename, ios::out | ios::trunc);

                // Write Header
                vOut << "Np" << "," << "Wall-time"  << endl;
                vOut << size << "," << wallTime << endl;
                vOut << "x" << "," << "y" << "," << "s" << "," << "v" << "," << "velU" << "," << "velV" << endl;

                // (x,y) coordinates of current point
                double x,y;

                for (int i = xShift_Start; i < (Nx - xShift_End); i++ ){
                    for (int j = yShift_Start; j < (Ny - yShift_End); j++ ){

                        y = subPos[0] + dy * (j - yShift_Start);
                        x = subPos[1] + dx * (i - xShift_Start);

                        vOut << x << "," << y << "," << s[j + i*Ny] << "," << v[j + i*Ny] << "," << velU[j + i*Ny] << "," << velV[j + i*Ny] << endl;
                    }
                }
                vOut.close();
            }
            else if (k == rank){
                // Open file and set as output only and append
                ofstream vOut(filename, ios::out | ios::app);

                // (x,y) coordinates of current point
                double x,y;

                for (int i = xShift_Start; i < (Nx - xShift_End); i++ ){
                    for (int j = yShift_Start; j < (Ny - yShift_End); j++ ){

                        y = subPos[0] + dy * (j - yShift_Start);
                        x = subPos[1] + dx * (i - xShift_Start);

                        vOut << x << "," << y << "," << s[j + i*Ny] << "," << v[j + i*Ny] << "," << velU[j + i*Ny] << "," << velV[j + i*Ny] << endl;
                    }
                }

                vOut.close();
            }
            // Ensure file is accessed by a single process at a time
            MPI_Barrier(MPIcomm);
        }

    }
