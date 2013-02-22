#ifndef HYDROGENICWAVEFUNCTION_H
#define HYDROGENICWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"

class HydrogenicWavefunction : public Wavefunction
{
public:
    HydrogenicWavefunction();
    double wavefunction(const mat &r);
    double laplace(const mat &r);
    mat gradient(const mat &r);

private:
    double nFactor;




};

#endif // HYDROGENICWAVEFUNCTION_H

