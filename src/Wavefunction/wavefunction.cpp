#include "wavefunction.h"

Wavefunction::Wavefunction()
{
}


double Wavefunction::waveFunction(int nDimensions, int nParticles, const mat &r, const double &alpha,const double &beta)
{
    argument=0.0;

    for (int i=0; i<nParticles; i++) {
        rSingleParticle=0;
        for (int j=0; j<nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }

    double TrialWaveFunction = exp(-alpha* argument)*jastrowFactor(nDimensions,nParticles,r,beta);
    return TrialWaveFunction;

}


double Wavefunction::jastrowFactor(int nDimensions, int nParticles,const mat &r,const double &beta){
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
