#ifndef HYDROGENICWAVEFUNCTION_H
#define HYDROGENICWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"

class HydrogenicWavefunction : public Wavefunction
{
public:
    HydrogenicWavefunction(Config *cfg);
    double waveFunction(int nParticles, const mat &r);
    double laplace(int nParticles, const mat &r, Config *cfg);
    Config* cfg;

private:
    int charge;
    double nFactor;



};

#endif // HYDROGENICWAVEFUNCTION_H

