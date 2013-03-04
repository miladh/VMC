#include "jastrowwavefunction.h"

JastrowWavefunction::JastrowWavefunction()
{
}

/************************************************************
Name:               JastrowWavefunction
Description:        jastrow wavefunction wavefunction
*/

double JastrowWavefunction::wavefunction(const mat &r)
{
    TrialWavefunction=1.0;

    for(uint i=0; i< r.n_rows; i++){  //CHECK THIS LATER!!!!
        for (uint qNum = 0; qNum < r.n_rows/2; qNum++){
            TrialWavefunction *= orbitals->orbitalEvaluate(r,qNum,i);
        }
    }

    TrialWavefunction *=jas.evaluateJastrow(r);

    return TrialWavefunction;

}

/************************************************************
Name:          Gradient
Description:
*/
mat JastrowWavefunction::gradient(const mat &r){
    dwavefunction=zeros<mat>(r.n_rows,r.n_cols);
    dHydrogenic=zeros<mat>(r.n_rows,r.n_cols);
    dJastrow=zeros<mat>(r.n_rows,r.n_cols);

    if(useAnalyticGradient){
        for (uint i = 0; i < r.n_rows; i++){
            for (uint qNum = 0; qNum < r.n_rows/2; qNum++){
                dHydrogenic.row(i)=orbitals->gradientOrbitalEvaluate(r,qNum,i);
            }
            dJastrow.row(i)=jas.gradientJastrowEvaluate(r,i);
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
double JastrowWavefunction::laplace(const mat &r){

    if(useAnalyticLaplace){
        ddwavefunction = 0;
        for (uint i = 0; i < r.n_rows; i++){
            for (uint qNum = 0; qNum < r.n_rows/2; qNum++){
                ddwavefunction+=orbitals->laplaceOrbitalEvaluate(r,qNum,i)+
                2*dot(jas.gradientJastrowEvaluate(r,i),orbitals->gradientOrbitalEvaluate(r,qNum,i));
            }
        }

        ddwavefunction+= jas.laplaceJastrowEvaluate(r);

        return ddwavefunction;
    }
    else{
        return laplaceNumerical(r);
    }

}



