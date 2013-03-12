#ifndef MCIP_H
#define MCIP_H

#include "src/Solver/solver.h"

class MCIS : public Solver
{
public:
    MCIS(const uint &nParticles, const uint &nDimensions,Hamiltonian *hamiltonian, Wavefunction *TrialWavefunction);
    void solve(int nCycles, long idum);


private:
    void MetropolisAlgoIS();
    mat getQuantumForce(const mat &r);
    mat qForce, qForceOld, qForceNew;

    long idum;
    double timeStep,nCycles;
    double energySum,energySquaredSum, deltaE;
    double h,D,GreensFunction;
};

#endif // MCIP_H
