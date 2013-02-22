#ifndef MCBF_H
#define MCBF_H

#include "src/Solver/solver.h"


class MCBF : public Solver
{
public:
    MCBF(Hamiltonian *hamiltonian, Wavefunction *TrialWaveFunction);
    void solve(int nCycles, long idum);


private:
    void MetropolisAlgoBF(int nCycles, double stepLength);
    double optimalStepLength();

    Hamiltonian *hamiltonian;
    Wavefunction* TrialWaveFunction;

    mat rOld,rNew;

    long idum;
    double stepLength,stepMinMax,stepMin;
    double waveFunctionOld,waveFunctionNew;
    double energySum,energySquaredSum, deltaE;

};

#endif // MCBF_H

