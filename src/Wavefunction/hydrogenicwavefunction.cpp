#include "hydrogenicwavefunction.h"

HydrogenicWavefunction::HydrogenicWavefunction()
{
}

/************************************************************
Name:               HydrogenicWavefunction
Description:        hydrogen like wavefunction
*/

double HydrogenicWavefunction::wavefunction(const mat &r)
{

    orbitals->k=charge;
    nFactor=4*pow((charge/sqrt(4*acos(-1))),3);
    TrialWavefunction = orbitals->orbitalEvaluate(r,0,0)*orbitals->orbitalEvaluate(r,0,1);

    return nFactor*TrialWavefunction;

}



/************************************************************
Name:          Gradient
Description:
*/
mat HydrogenicWavefunction::gradient(const mat &r){
    orbitals->k=charge;
    nFactor=4*pow((charge/sqrt(4*acos(-1))),3);

    if(useAnalyticGradient){
     for (uint i = 0; i < r.n_rows; i++){
            dwavefunction.row(i)=orbitals->gradientOrbitalEvaluate(r,0,i);
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
double HydrogenicWavefunction::laplace(const mat &r){
    orbitals->k=charge;
    nFactor=4*pow((charge/sqrt(4*acos(-1))),3);

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




