#include "hlikewavefunction.h"
#include"src/Jastrow/nojastrow.h"
#include"src/Jastrow/padejastrow.h"
#include"src/Orbitals/hydrogenic.h"


HLikeWavefunction::HLikeWavefunction(const uint &nParticles):
    Wavefunction(nParticles),
    dHydrogenic(zeros(nParticles,nParticles)),
    dJastrow(zeros(nParticles,nParticles))
{
}


/************************************************************
Name:
Description:
*/

double HLikeWavefunction::wavefunction(const mat &r)
{
    TrialWavefunction =slater->initializeSD(r)*jas->evaluateJastrow(r);

    return TrialWavefunction;

}

/************************************************************
Name:          Gradient
Description:
*/
mat HLikeWavefunction::gradient(const mat &r){

    if(useAnalyticGradient){
        for (uint i = 0; i < r.n_rows; i++){
            for (uint qNum = 0; qNum < r.n_rows/2; qNum++){
                dHydrogenic.row(i)=orbitals->gradientOrbitalEvaluate(r,qNum,i);
            }
            dJastrow.row(i)=jas->gradientJastrowEvaluate(r,i);
        }

        dwavefunction=dHydrogenic+dJastrow;
        return dwavefunction;
    }
    else{
        return gradientNumerical(r);
    }

}


/************************************************************
Name:          laplace
Description:
*/
double HLikeWavefunction::laplace(const mat &r){

    if(useAnalyticLaplace){
        ddwavefunction = 0;
        for (uint i = 0; i < r.n_rows; i++){
            for (uint qNum = 0; qNum < r.n_rows/2; qNum++){
                ddwavefunction+=orbitals->laplaceOrbitalEvaluate(r,qNum,i)+
                2*dot(jas->gradientJastrowEvaluate(r,i),orbitals->gradientOrbitalEvaluate(r,qNum,i));
            }
        }

        ddwavefunction+= jas->laplaceJastrowEvaluate(r);

        return ddwavefunction;
    }
    else{
        return laplaceNumerical(r);
    }

}



