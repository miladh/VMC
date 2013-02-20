#ifndef MCIP_H
#define MCIP_H

#include "src/Solver/solver.h"

class MCIS : public Solver
{
public:
    MCIS(Hamiltonian *hamiltonian, Wavefunction *TrialWaveFunction);
    void solve(int nCycles, long idum);
    void loadConfiguration(Config *cfg);

private:


    int nDimensions,nParticles;
    int nCycles,thermalization;
    long idum;
    double timeStep;

    Hamiltonian *hamiltonian;
    Wavefunction* TrialWaveFunction;

    mat rOld;
    mat rNew;

    double waveFunctionOld,waveFunctionNew;
    double energySum,energySquaredSum, deltaE;

    void MetropolisAlgoIS(int nCycles, long idum);
    mat QuantumForce(const mat &r);


    mat qForceOld;
    mat qForceNew;
    mat qForce;
    double h,D,GreensFunction;

};

#endif // MCIP_H
