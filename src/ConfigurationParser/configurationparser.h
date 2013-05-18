#ifndef CONFIGURATIONPARSER_H
#define CONFIGURATIONPARSER_H

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

class ConfigurationParser
{
public:
    ConfigurationParser(Config *cfg, const int &myRank, const int &nProcess);

    void runSolver();
    double getEnergy();
    double getEnergySquared();
    double getVariance();
    double getSigma();
    double getAcceptanceRate();
    vec getVariationalDerivate();

    double alpha, beta;

    void setVariationalParameters(vector<double> paramters);


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
    void setup();
    void setWavefunction();
    void setHamiltonian();
    void setObservables();
    void setSolver();


};

#endif // CONFIGURATIONPARSER_H


