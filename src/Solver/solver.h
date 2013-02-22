#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>

#include "src/Wavefunction/wavefunction.h"
#include "src/Hamiltonian/hamiltonian.h"
#include "src/Potential/potential.h"
#include "src/Kinetic/kinetic.h"

using namespace arma;
using namespace std;
using namespace libconfig;


class Solver
{
public:
    Solver();
    virtual void solve(int nCycles, long idum) = 0;
    void loadConfiguration(Config *cfg);

    double energySquared;
    double energy;
    double acceptedSteps;


protected:
    int nDimensions,nParticles;
    double thermalization,nPreCycles;
    double minStepLength, maxStepLength, tolerance;
};

#endif // VMCSOLVER_H
