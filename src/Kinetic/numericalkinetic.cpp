#include "numericalkinetic.h"
#include "src/Wavefunction/jastrowwavefunction.h"
#include "src/Wavefunction/basicwavefunction.h"

NumericalKinetic::NumericalKinetic():
    h(0.001),
    h2(1000000)
{
}

double NumericalKinetic::evaluate(int nParticles,const mat &r){

    hVec=h*ones<rowvec>(r.n_cols);
    rPlus = zeros<mat>(nParticles, r.n_cols);
    rMinus = zeros<mat>(nParticles, r.n_cols);

    rPlus = rMinus = r;

    waveFunctionMinus = 0;
    waveFunctionPlus = 0;

    waveFunctionCurrent = wf->waveFunction(nParticles,r);

    KineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {                
            rPlus.row(i) += hVec;
            rMinus.row(i) -= hVec;
            waveFunctionMinus = wf->waveFunction(nParticles,rMinus);
            waveFunctionPlus = wf->waveFunction(nParticles,rPlus);
            KineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus.row(i) = r.row(i);
            rMinus.row(i)= r.row(i);
        }

    KineticEnergy = 0.5 * h2 * KineticEnergy / waveFunctionCurrent;

    return KineticEnergy;
}
