#ifndef MCBF_H
#define MCBF_H

#include "src/Solver/solver.h"


class MCBF : public Solver
{
public:
    MCBF();
    void solve(int nDimensions,int nParticles, Hamiltonian* hamiltonian, Wavefunction* TrialWaveFunction, int nCycles, long idum);

private:
    double stepLength;
    mat rOld;
    mat rNew;

};

#endif // MCBF_H

