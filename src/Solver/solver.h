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
    virtual void solve() = 0;
    virtual void loadConfiguration(Config *cfg)=0;

    double energySquared;
    double energy;

};

#endif // VMCSOLVER_H
