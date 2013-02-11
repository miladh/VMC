#ifndef HYDROGENICWAVEFUNCTION_H
#define HYDROGENICWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"

class HydrogenicWavefunction : public Wavefunction
{
public:
    HydrogenicWavefunction();
    double waveFunction(int nParticles, const mat &r);

private:
    int charge;
    double nFactor;
    double rSingleParticle;
    double argument;


};

#endif // HYDROGENICWAVEFUNCTION_H

