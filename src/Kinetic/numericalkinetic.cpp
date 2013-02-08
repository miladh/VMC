#include "numericalkinetic.h"
#include "src/Wavefunction/jastrowwavefunction.h"
#include "src/Wavefunction/basicwavefunction.h"

NumericalKinetic::NumericalKinetic():
    h(0.001),
    h2(1000000)
{
}

double NumericalKinetic::evaluate(int nDimensions, int nParticles,const mat &r){

    rPlus = zeros<mat>(nParticles, nDimensions);
    rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    waveFunctionMinus = 0;
    waveFunctionPlus = 0;

    waveFunctionCurrent = wf->waveFunction(nDimensions,nParticles,r);



    KineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = wf->waveFunction(nDimensions,nParticles,rMinus);
            waveFunctionPlus = wf->waveFunction(nDimensions,nParticles,rPlus);
            KineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }


    KineticEnergy = 0.5 * h2 * KineticEnergy / waveFunctionCurrent;

    return KineticEnergy;
}
