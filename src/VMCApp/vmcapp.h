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
    VMCApp(Config *cfg, const int &myRank, const int &nProcess);

    void runVMCApp();
    double getEnergy();
    double getEnergySquared();
    double getVariance();
    double getSigma();
    double getAcceptanceRate();
    vec getVariationalDerivate();

    double alpha, beta;


private:
    Config* cfg;
    Orbitals* orbitals;
    Jastrow* jastrow;
    Wavefunction* trialWavefunction;
    Solver* solver;
    Hamiltonian *hamiltonian;
    Observables* observables;

    int nParticles,nDimensions,charge;
    int nProcess, myRank;
    int systemType, wavefunctionType,solverType,InteractionType;
    int minimizationIsEnable, blockingIsEnable ;
    long idum;
    double nCycles;

    double R;
    double totEnergy,totEnergySquared,averageDistance;
    double Variance, Acceptance,Sigma;
    vec totVariationalDerivate,totEnergyVarDerivate;
    double tmp;
    vec tmpVec;





    void loadAndSetConfiguration();
    void setWavefunction();
    void setHamiltonian();
    void setObservables();
    void setAndRunSolver();
    void messagePassing();


};

#endif // VMCAPP_H
