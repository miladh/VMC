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

    TrialWavefunction = orbitals->orbitalEvaluate(r,0,0)*orbitals->orbitalEvaluate(r,0,1);
    TrialWavefunction *=jas.JastrowExponential(r.n_rows,r);

    return TrialWavefunction;

}

/************************************************************
Name:          Gradient
Description:
*/
mat JastrowWavefunction::gradient(const mat &r){

    if(useAnalyticGradient){
        for (uint i = 0; i < r.n_rows; i++){
            dHydrogenic.row(i)=orbitals->GradientOrbitalEvaluate(r,0,i);
        }

        return dHydrogenic;
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
                ddHydrogenic += orbitals->LaplaceOrbitalEvaluate(r,qNum,i); //*SlaterInv(j, i)
            }
        }
        ddwavefunction= ddHydrogenic+jas.LaplaceJastrowEvaluate(r);

        return ddwavefunction;
    }
    else{
        return laplaceNumerical(r);
    }

}



