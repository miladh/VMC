#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

#include "src/Wavefunction/wavefunction.h"
#include "src/Hamiltonian/hamiltonian.h"
#include "src/Potential/potential.h"
#include "src/Kinetic/kinetic.h"
#include "src/includes/Defines.h"

using namespace arma;
using namespace std;
using namespace libconfig;


class Solver
{
public:
    Solver(const uint &nParticles, const uint &nDimensions, Hamiltonian *hamiltonian, Wavefunction *TrialWavefunction);
    virtual void solve(int nCycles, long idum) = 0;
    void loadConfiguration(Config *cfg);

    double energySquared;
    double energy;
    double acceptedSteps;


protected:
    int nParticles,nDimensions;
    double thermalization,nPreCycles;
    double minStepLength, maxStepLength, tolerance;
    double R;
    mat rOld, rNew;
    Hamiltonian *hamiltonian;
    Wavefunction *TrialWavefunction;

    int myRank,nProcess;
    mat energyVector;

};

#endif // VMCSOLVER_H
