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
    void writeToFile(ofstream myfile);
    double getEnergy();
    double getEnergySquared();
    double getVariance();
    double getSigma();
    double getAcceptanceRate();

    Wavefunction* setWaveFunction();
    Solver* setSolverMethod();


    Config *cfg;
    Solver* solver;
    Wavefunction *TrialWaveFunction;
    Potential *potential;
    Kinetic *kinetic;
    Hamiltonian *hamiltonian;
    double alpha, beta;

private:
    int nDimensions,nParticles,nCycles;
    int nProcess, myRank;
    long idum;

    double totEnergy,totEnergySquared;
    double Variance, Acceptance,Sigma;
    double tmp;

    ofstream myfile;

};

#endif // VMCAPP_H
