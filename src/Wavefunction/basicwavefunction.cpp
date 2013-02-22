#include "basicwavefunction.h"

BasicWavefunction::BasicWavefunction()
{
}


/************************************************************
Name:               BasicWaveFunction
Description:        simple wavefunction
*/

double BasicWavefunction::wavefunction(const mat &r)
{
    TrialWavefunction = orbitals->orbitalEvaluate(r,0,0)*orbitals->orbitalEvaluate(r,0,1);

    return TrialWavefunction;
}


/************************************************************
Name:          Gradient
Description:
*/
mat BasicWavefunction::gradient(const mat &r){

    if(useAnalyticGradient){
        for (uint i = 0; i < r.n_rows; i++){
            dwavefunction.row(i)=orbitals->GradientOrbitalEvaluate(r,0,i);
        }
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
double BasicWavefunction::laplace(const mat &r){

    if(useAnalyticLaplace){
        ddwavefunction = 0;
        for (uint i = 0; i < r.n_rows; i++){
            for (uint qNum = 0; qNum < r.n_rows/2; qNum++){
                ddwavefunction += orbitals->LaplaceOrbitalEvaluate(r,qNum,i); //*SlaterInv(j, i)
            }
        }
        return ddwavefunction;
    }
    else{
        return laplaceNumerical(r);
    }

}



