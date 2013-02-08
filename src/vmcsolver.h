#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "../src/Wavefunction/wavefunction.h"

using namespace arma;

class VMCSolver
{
public:
    VMCSolver();

    void runMonteCarloIntegration();

private:
    Wavefunction TrialWaveFunction;

    double localEnergy(const mat &r);

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;

    double h;
    double h2;

    long idum;

    double alpha;
    double beta;

    int nCycles;

    mat rOld;
    mat rNew;
};

#endif // VMCSOLVER_H
