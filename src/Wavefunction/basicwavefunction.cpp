#include "basicwavefunction.h"

BasicWaveFunction::BasicWaveFunction()
{
}



double BasicWaveFunction::waveFunction(int nDimensions, int nParticles, const mat &r)
{
    argument=0.0;

    for (int i=0; i<nParticles; i++) {
        rSingleParticle=0;
        for (int j=0; j<nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }

    double TrialWaveFunction = exp(-alpha* argument);
    return TrialWaveFunction;

}
