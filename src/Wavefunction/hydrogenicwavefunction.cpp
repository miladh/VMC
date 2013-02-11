#include "hydrogenicwavefunction.h"

HydrogenicWavefunction::HydrogenicWavefunction():
    charge(2),
    nFactor(4*pow((charge/sqrt(4*acos(-1))),3))
{
}

double HydrogenicWavefunction::waveFunction(int nParticles, const mat &r)
{
    argument=0.0;

    for (int i=0; i<nParticles; i++) {
        rSingleParticle = norm(r.row(i),2);
        argument += rSingleParticle;
    }

    TrialWaveFunction = exp(-charge* argument);
    return TrialWaveFunction*nFactor;

}
