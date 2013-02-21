#include "jastrowwavefunction.h"

JastrowWavefunction::JastrowWavefunction()
{
}

/************************************************************
Name:               JastrowWavefunction
Description:        jastrow wavefunction wavefunction
*/

double JastrowWavefunction::waveFunction(int nParticles, const mat &r)
{

    TrialWaveFunction = orbitals->orbitalEvaluate(r,0,0)*orbitals->orbitalEvaluate(r,0,1);
    TrialWaveFunction *=jas.JastrowExponential(nParticles,r);


    return TrialWaveFunction;

}


/************************************************************
Name:          laplace
Description:
*/
double JastrowWavefunction::laplace(int nParticles, const mat &r, Config* cfg){

    analytic= cfg->lookup("AppSettings.useAnalyticLaplace");

    if(analytic){

        ddHydrogenic = 0;
        for (int i = 0; i < nParticles; i++){
            for (int qNum = 0; qNum < nParticles/2; qNum++){
                ddHydrogenic += orbitals->LaplaceOrbitalEvaluate(r,qNum,i); //*SlaterInv(j, i)
            }
        }

        ddwaveFunction= ddHydrogenic+jas.LaplaceJastrowEvaluate(r);

        return ddwaveFunction;

    }

    else{
        return laplaceNumerical(nParticles,r,cfg);
    }

}





