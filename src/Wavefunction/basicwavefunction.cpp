#include "basicwavefunction.h"

BasicWaveFunction::BasicWaveFunction()
{
}



double BasicWaveFunction::waveFunction( int nParticles, const mat &r)
{
    argument=0.0;

    for (int i=0; i<nParticles; i++) {
        rSingleParticle = norm(r.row(i),2);
        argument += rSingleParticle;
    }

    TrialWaveFunction = exp(-alpha* argument);
    return TrialWaveFunction;

}
