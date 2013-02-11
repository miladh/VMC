#ifndef MCBF_H
#define MCBF_H

#include "src/Solver/solver.h"


class MCBF : public Solver
{
public:
    MCBF(Hamiltonian *hamiltonian, Wavefunction *TrialWaveFunction);
    void solve();
    void loadConfiguration(Config *cfg);

private:
    int nDimensions,nParticles;
    int nCycles,nPreCycles;
    long idum;
    double stepLength,acceptedSteps;

    Hamiltonian *hamiltonian;
    Wavefunction* TrialWaveFunction;

    mat rOld;
    mat rNew;

    double waveFunctionOld,waveFunctionNew;
    double energySum,energySquaredSum, deltaE;
    double minStepLength, maxStepLength, tolerance;

    void MetropolisAlgo(int nCycles, double stepLength);
    double optimalStepLength();
};

#endif // MCBF_H

