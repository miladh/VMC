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
    jas.setaValues(r.n_rows);
    TrialWavefunction = orbitals->orbitalEvaluate(r,0,0)*orbitals->orbitalEvaluate(r,0,1);
    TrialWavefunction *=jas.evaluateJastrow(r);

    return TrialWavefunction;

}

/************************************************************
Name:          Gradient
Description:
*/
mat JastrowWavefunction::gradient(const mat &r){

    if(useAnalyticGradient){
        for (uint i = 0; i < r.n_rows; i++){
            dHydrogenic.row(i)=orbitals->gradientOrbitalEvaluate(r,0,i);
            dJastrow.row(i)=jas.gradientJastrowEvaluate(r,i);
        }

        return (dHydrogenic+dJastrow);
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
        ddHydrogenic = 0;
        for (uint i = 0; i < r.n_rows; i++){
            for (uint qNum = 0; qNum < r.n_rows/2; qNum++){
                ddHydrogenic += orbitals->laplaceOrbitalEvaluate(r,qNum,i); //*SlaterInv(j, i)
            }
        }
        ddwavefunction= ddHydrogenic+jas.laplaceJastrowEvaluate(r)+
                2*dot(jas.gradientJastrowEvaluate(r,0),orbitals->gradientOrbitalEvaluate(r,0,0))+
                2*dot(jas.gradientJastrowEvaluate(r,1),orbitals->gradientOrbitalEvaluate(r,0,1));


        return ddwavefunction;
    }
    else{
        return laplaceNumerical(r);
    }

}



