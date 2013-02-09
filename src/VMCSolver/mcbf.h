#ifndef MCBF_H
#define MCBF_H

#include "src/VMCSolver/vmcsolver.h"


class MCBF : public VMCSolver
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

