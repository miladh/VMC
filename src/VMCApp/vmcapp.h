#ifndef VMCAPP_H
#define VMCAPP_H


#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include "src/Solver/solver.h"
#include "src/Wavefunction/wavefunction.h"
#include "src/Hamiltonian/hamiltonian.h"

using namespace arma;
using namespace std;
using namespace libconfig;

class VMCApp
{
public:
    VMCApp();

    Solver* solver;
    Wavefunction *TrialWaveFunction;
    Potential *potential;
    Kinetic *kinetic;
    Hamiltonian *hamiltonian;

    int nDimensions,nParticles,nCycles;
    double alpha, beta;
    long idum;

    double energySquared;
    double energy;

    mat rOld;
    mat rNew;

    void runVMCApp();
};

#endif // VMCAPP_H
