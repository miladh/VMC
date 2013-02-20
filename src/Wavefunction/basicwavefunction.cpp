#include "basicwavefunction.h"

BasicWavefunction::BasicWavefunction()
{
}


/************************************************************
Name:               BasicWaveFunction
Description:        simple wavefunction
*/

double BasicWavefunction::waveFunction( int nParticles, const mat &r)
{
    argument=0.0;

    for (int i=0; i<nParticles; i++) {
        rSingleParticle = norm(r.row(i),2);
        argument += rSingleParticle;
    }

    TrialWaveFunction = exp(-alpha* argument);
    return TrialWaveFunction;

}


/************************************************************
Name:          laplace
Description:
*/
double BasicWavefunction::laplace(int nParticles, const mat &r, Config* cfg){
    analytic= cfg->lookup("AppSettings.useAnalyticLaplace");

    if(analytic){
    r1 = norm(r.row(0), 2);
    r2 = norm(r.row(1), 2);
    rij = norm(r.row(0) - r.row(1), 2);

    ddwaveFunction = 2*alpha*(alpha-1/r1-1/r2);

    return ddwaveFunction;
    }else{
        return laplaceNumerical(nParticles,r,cfg);
    }

}


