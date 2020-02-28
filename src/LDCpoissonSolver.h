#pragma once

#include <string>
#include <cstring>
using namespace std;

class LDCpoissonSolver
{
    public:
    // Constructor
    LDCpoissonSolver();
    ~LDCpoissonSolver();

    void Initialize();

    private:
    // Linear solver variables
    double* A = nullptr;

};