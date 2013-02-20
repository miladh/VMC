#ifndef MCBF_H
#define MCBF_H

#include "src/Solver/solver.h"


class MCBF : public Solver
{
public:
    MCBF(Hamiltonian *hamiltonian, Wavefunction *TrialWaveFunction);
    void solve(int nCycles, long idum);
    void loadConfiguration(Config *cfg);

private:
    int nDimensions,nParticles;
    int nCycles,nPreCycles,thermalization;
    long idum;
    double stepLength;

    Hamiltonian *hamiltonian;
    Wavefunction* TrialWaveFunction;

    mat rOld;
    mat rNew;

    double waveFunctionOld,waveFunctionNew;
    double energySum,energySquaredSum, deltaE;
    double minStepLength, maxStepLength, tolerance;

    void MetropolisAlgoBF(int nCycles, double stepLength, long idum);
    double optimalStepLength(long idum);
};

#endif // MCBF_H

