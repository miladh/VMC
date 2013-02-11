#include "jastrowwavefunction.h"

JastrowWavefunction::JastrowWavefunction()
{
}

double JastrowWavefunction::waveFunction(int nParticles, const mat &r)
{
    argument=0.0;

    for (int i=0; i<nParticles; i++) {
        rSingleParticle = norm(r.row(i),2);
        argument += rSingleParticle;
    }

    TrialWaveFunction = exp(-alpha* argument);

    TrialWaveFunction *=jastrowFactor(nParticles,r);
    return TrialWaveFunction;

}

double JastrowWavefunction::jastrowFactor(int nParticles, const mat &r){
    correlation=0;

    for (int i=0; i<nParticles-1; i++) {
        for (int j=i+1; j<nParticles; j++) {
            rij= norm( r.row(i)-r.row(j) ,2);

            correlation+=rij/(2+2*beta*rij);
        }
    }

    double Jastrow = exp(correlation);
    return Jastrow;

}
