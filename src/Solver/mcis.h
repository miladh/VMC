#ifndef MCIP_H
#define MCIP_H

#include <src/Solver/solver.h>

class MCIS : public Solver
{
public:
    MCIS(Config* cfg,Hamiltonian *hamiltonian, Wavefunction *TrialWavefunction, Observables *observables);
    void solve(int nCycles, long idum);


private:
    void MetropolisAlgoIS();
    mat getQuantumForce(const mat &r);
    mat qForce, qForceOld, qForceNew;

    long idum;
    double nCycles;
    double deltaE;
    double GreensFunction;
    int myRank;
};

#endif // MCIP_H
