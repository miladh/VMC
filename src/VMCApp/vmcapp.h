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
    VMCApp(Config* cfg);
    void runVMCApp(int nCycles, long idum);
    void writeToFile(ofstream myfile);

    Config *cfg;
    Solver* solver;
    Wavefunction *TrialWaveFunction;
    Potential *potential;
    Kinetic *kinetic;
    Hamiltonian *hamiltonian;


    int nDimensions,nParticles,nCycles;
    int nProcess, myRank;
    long idum;
    double alpha, beta;
    double energy,energySquared,totEnergy,totEnergySquared;
    double Variance, Acceptance,Sigma;
    double tmp;

    ofstream myfile;

private:
    mat rOld;
    mat rNew;

};

#endif // VMCAPP_H
