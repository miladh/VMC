#ifndef MCBF_H
#define MCBF_H

#include <src/Solver/solver.h>


class MCBF : public Solver
{
public:
    MCBF(const uint &nParticles, const uint &nDimensions,Hamiltonian *hamiltonian, Wavefunction *TrialWavefunction);
    void solve(int nCycles, long idum);


private:
    void MetropolisAlgoBF(int nCycles, double stepLength);
    double optimalStepLength();

    long idum;
    double stepLength,stepMinMax,stepMin;
    double  deltaE;


};

#endif // MCBF_H

