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
    orbitals->k=charge;
    TrialWaveFunction = orbitals->orbitalEvaluate(r,0,0)*orbitals->orbitalEvaluate(r,0,1);

    return nFactor*TrialWaveFunction;

}


/************************************************************
Name:          laplace
Description:
*/
double HydrogenicWavefunction::laplace(int nParticles, const mat &r, Config* cfg){
    analytic= cfg->lookup("AppSettings.useAnalyticLaplace");
    orbitals->k=charge;

    if(analytic){

        ddwaveFunction = 0;
        for (int i = 0; i < nParticles; i++){
            for (int qNum = 0; qNum < nParticles/2; qNum++){
                ddwaveFunction += orbitals->LaplaceOrbitalEvaluate(r,qNum,i); //*SlaterInv(j, i)
            }
        }

        return ddwaveFunction;
    }

    else{
        return laplaceNumerical(nParticles,r,cfg);
    }

}



