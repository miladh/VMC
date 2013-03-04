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
    TrialWavefunction=slaterDet.evaluateSlater(r);

    return TrialWavefunction;
}


/************************************************************
Name:          Gradient
Description:
*/
mat BasicWavefunction::gradient(const mat &r){
    dwavefunction=zeros<mat>(r.n_rows,r.n_cols);

    if(useAnalyticGradient){
        for (uint i = 0; i < r.n_rows; i++){
            for (uint qNum = 0; qNum < r.n_rows/2; qNum++){
                dwavefunction.row(i)=orbitals->gradientOrbitalEvaluate(r,qNum,i);
            }
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
                ddwavefunction += orbitals->laplaceOrbitalEvaluate(r,qNum,i); //*SlaterInv(j, i)
            }
        }
        return ddwavefunction;
    }
    else{
        return laplaceNumerical(r);
    }

}



