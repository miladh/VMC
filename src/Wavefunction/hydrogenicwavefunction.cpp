#include "hydrogenicwavefunction.h"

HydrogenicWavefunction::HydrogenicWavefunction(Config* cfg):
    charge(cfg->lookup("PotentialSettings.charge")),
    nFactor(4*pow((charge/sqrt(4*acos(-1))),3))
{
    this->cfg=cfg;
}



/************************************************************
Name:               HydrogenicWavefunction
Description:        hydrogen like wavefunction
*/

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


/************************************************************
Name:          laplace
Description:
*/
double HydrogenicWavefunction::laplace(int nParticles, const mat &r, Config* cfg){
    analytic= cfg->lookup("AppSettings.useAnalyticLaplace");

    if(analytic){
    r1 = norm(r.row(0), 2);
    r2 = norm(r.row(1), 2);
    rij = norm(r.row(0) - r.row(1), 2);

    ddwaveFunction = 2*charge*(charge-1/r1-1/r2);

    return ddwaveFunction;
    }else{
        return laplaceNumerical(nParticles,r,cfg);
    }
}



