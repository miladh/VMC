#ifndef JASTROWWAVEFUNCTION_H
#define JASTROWWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"
#include "src/Jastrow/jastrow.h"

class JastrowWavefunction : public Wavefunction
{
public:
    JastrowWavefunction();

    double jastrowFactor(const mat &r);
    double wavefunction(const mat &r);
    double laplace(const mat &r);
    mat gradient(const mat &r);

private:
    mat dHydrogenic,dJastrow;

};

#endif // JASTROWWAVEFUNCTION_H
