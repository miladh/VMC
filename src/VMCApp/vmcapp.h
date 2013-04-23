#ifndef VMCAPP_H
#define VMCAPP_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include <src/Solver/solver.h>
#include <src/Wavefunction/wavefunction.h>
#include <src/Hamiltonian/hamiltonian.h>
#include <src/electronInteraction/electroninteraction.h>
#include <src/Observables/observables.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
#pragma GCC diagnostic warning "-Wunused-parameter"

using namespace arma;
using namespace std;
using namespace libconfig;

class VMCApp
{
public:
    VMCApp(const int &myRank, const int &nProcess);

    void loadConfiguration(Config *cfg);
    void runVMCApp(int nCycles, long idum);
    void messagePassing();
    double getEnergy();
    double getEnergySquared();
    double getVariance();
    double getSigma();
    double getAcceptanceRate();
    vec getVariationalDerivate();

    double alpha, beta;


private:
    int nParticles,nDimensions,charge;
    int nProcess, myRank;
    int WavefunctionType,solverType,InteractionType;
    int minimizationIsEnable, blockingIsEnable ;
    long idum;

    double totEnergy,totEnergySquared;
    double Variance, Acceptance,Sigma;
    vec totVariationalDerivate,totEnergyVarDerivate;
    double tmp;
    vec tmpVec;

    Config *cfg;
    Solver* solver;
    Wavefunction *trialWavefunction;
    Potential *potential;
    Kinetic *kinetic;
    ElectronInteraction *electonInteraction;
    Hamiltonian *hamiltonian;
    Observables* observables;

    Wavefunction* setWavefunction();
    Solver* setSolverMethod();
    ElectronInteraction* setInteraction();

};

#endif // VMCAPP_H
