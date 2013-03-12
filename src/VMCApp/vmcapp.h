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
    VMCApp(Config* cfg, const int &myRank, const int &nProcess);

    void runVMCApp(int nCycles, long idum);
    double getEnergy();
    double getEnergySquared();
    double getVariance();
    double getSigma();
    double getAcceptanceRate();

    double alpha, beta;


private:
    int nParticles,charge,nDimensions;
    int nProcess, myRank;
    int WavefunctionType,solverType;
    long idum;

    double totEnergy,totEnergySquared;
    double Variance, Acceptance,Sigma;
    double tmp;

    Config *cfg;
    Solver* solver;
    Wavefunction *TrialWavefunction;
    Potential *potential;
    Kinetic *kinetic;
    Hamiltonian *hamiltonian;

    Wavefunction* setWavefunction();
    Solver* setSolverMethod();

};

#endif // VMCAPP_H
