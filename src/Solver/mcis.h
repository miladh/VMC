#ifndef MCIP_H
#define MCIP_H

#include "src/Solver/solver.h"

class MCIS : public Solver
{
public:
    MCIS(Hamiltonian *hamiltonian, Wavefunction *TrialWavefunction);
    void solve(int nCycles, long idum);


private:
    void MetropolisAlgoIS();
    mat getQuantumForce(const mat &r);

    Hamiltonian *hamiltonian;
    Wavefunction* TrialWavefunction;

    mat qForce, qForceOld, qForceNew;
    mat rOld, rNew;

    long idum;
    double timeStep,nCycles;

    double wavefunctionOld,wavefunctionNew;
    double energySum,energySquaredSum, deltaE;
    double h,D,GreensFunction;
};

#endif // MCIP_H
