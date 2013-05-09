#ifndef MCIP_H
#define MCIP_H

#include <src/Solver/solver.h>

class MCIS : public Solver
{
public:
    MCIS(Hamiltonian *hamiltonian, Wavefunction *TrialWavefunction, Observables *observables);
    void solve(int nCycles, long idum);


private:
    void MetropolisAlgoIS();
    mat getQuantumForce(const mat &r);
    mat qForce, qForceOld, qForceNew;

    long idum;
    double timeStep,nCycles;
    double deltaE;
    double D,GreensFunction;
};

#endif // MCIP_H
