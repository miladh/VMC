#include "wavefunction.h"

Wavefunction::Wavefunction()
{
 orbitals = new Hydrogenic;
}


double Wavefunction::laplaceNumerical(int nParticles, const mat &r, Config *cfg){
    h=cfg->lookup("NumericalKineticSettings.h");
    h2=cfg->lookup("NumericalKineticSettings.h2");

    hVec=h*ones<rowvec>(r.n_cols);
    rPlus = zeros<mat>(nParticles, r.n_cols);
    rMinus = zeros<mat>(nParticles, r.n_cols);

    rPlus = rMinus = r;

    waveFunctionMinus = 0;
    waveFunctionPlus = 0;

    waveFunctionCurrent = waveFunction(nParticles,r);

    ddwaveFunction = 0;
    for(int i = 0; i < nParticles; i++) {
            rPlus.row(i) += hVec;
            rMinus.row(i) -= hVec;
            waveFunctionMinus = waveFunction(nParticles,rMinus);
            waveFunctionPlus = waveFunction(nParticles,rPlus);
            ddwaveFunction += (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus.row(i) = r.row(i);
            rMinus.row(i)= r.row(i);
        }

    ddwaveFunction =h2 * ddwaveFunction / waveFunctionCurrent;

    return ddwaveFunction;
}
