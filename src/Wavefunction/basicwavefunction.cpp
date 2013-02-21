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

    TrialWaveFunction = orbitals->orbitalEvaluate(r,0,0)*orbitals->orbitalEvaluate(r,0,1);

    return TrialWaveFunction;

}


/************************************************************
Name:          laplace
Description:
*/
double BasicWavefunction::laplace(int nParticles, const mat &r, Config* cfg){

    analytic= cfg->lookup("AppSettings.useAnalyticLaplace");

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


