#ifndef MCBF_H
#define MCBF_H

#include "src/Solver/solver.h"


class MCBF : public Solver
{
public:
    MCBF(Hamiltonian *hamiltonian, Wavefunction *TrialWaveFunction);
    void solve();

private:
    int nDimensions;
    int nParticles;
    Hamiltonian *hamiltonian;
    Wavefunction* TrialWaveFunction;
    int nCycles,nPreCycles;
    long idum;


    double stepLength,acceptedSteps;
    mat rOld;
    mat rNew;

    double waveFunctionOld ;
    double waveFunctionNew;

    double energySum ;
    double energySquaredSum;

    double deltaE;


    void MetropolisAlgo(int nCycles, double stepLength);
    double optimalStepLength();
    double difference(double stepLength);

};

#endif // MCBF_H

