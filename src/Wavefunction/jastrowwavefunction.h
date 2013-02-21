#ifndef JASTROWWAVEFUNCTION_H
#define JASTROWWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"
#include "src/Jastrow/jastrow.h"

class JastrowWavefunction : public Wavefunction
{
public:
    JastrowWavefunction();


    double jastrowFactor(int nParticles, const mat &r);
    double waveFunction(int nParticles, const mat &r);
    double laplace(int nParticles, const mat &r, Config *cfg);

private:
    double ddHydrogenic;





};

#endif // JASTROWWAVEFUNCTION_H
