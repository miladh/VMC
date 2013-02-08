#include "jastrowwavefunction.h"

JastrowWaveFunction::JastrowWaveFunction()
{
}

double JastrowWaveFunction::waveFunction(int nDimensions, int nParticles, const mat &r)
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

    TrialWaveFunction *=jastrowFactor(nDimensions,nParticles,r);
    return TrialWaveFunction;

}

double JastrowWaveFunction::jastrowFactor(int nDimensions, int nParticles, const mat &r){
    correlation=0;

    for (int i=0; i<nParticles-1; i++) {
        for (int j=i+1; j<nParticles; j++) {
            rij= 0.0;

            for (int k=0; k <nDimensions ; k++ ){
                rij += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
            }
            correlation+=sqrt(rij)/(2+2*beta*sqrt(rij));
        }
    }

    double Jastrow = exp(correlation);
    return Jastrow;

}
