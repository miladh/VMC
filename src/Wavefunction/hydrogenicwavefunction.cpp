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
