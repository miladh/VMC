#include "wavefunction.h"


Wavefunction::Wavefunction(const uint &nParticles) :
    nParticles(nParticles),
    dwavefunction(zeros<mat>(nParticles,3)),
    dvariational(zeros<vec>(2)),
    hGrad(1e-5),
    h(1e-5),
    h2(1e10),
    rPlus(zeros<mat>(nParticles, 3)),
    rMinus(zeros<mat>(nParticles, 3)),
    slater(new Slater(nParticles))
{
}

/************************************************************
Name:               laplaceNumerical
Description:
*/
double Wavefunction::laplaceNumerical(const mat &r){

    rPlus = rMinus = r;
    wavefunctionCurrent = wavefunction(r);
    ddwavefunction = 0;

    for(uint i = 0; i <r.n_rows; i++) {
        for(uint j=0; j <r.n_cols ; j++){
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            wavefunctionMinus = wavefunction(rMinus);
            wavefunctionPlus = wavefunction(rPlus);
            ddwavefunction += (wavefunctionMinus + wavefunctionPlus - 2 * wavefunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j)= r(i,j);
        }
    }

    ddwavefunction =h2 * ddwavefunction / wavefunctionCurrent;

    return ddwavefunction;
}



/************************************************************
Name:               gradientNumerical
Description:
*/
mat Wavefunction::gradientNumerical(const mat &r){

    rPlus = rMinus = r;

    wavefunctionCurrent = wavefunction(r);
    for(uint i =0; i<r.n_rows ; i++){
        for(uint j =0; j<r.n_cols ; j++){
            rPlus(i,j)+=hGrad;
            rMinus(i,j)-=hGrad;
            wavefunctionMinus = wavefunction(rMinus);
            wavefunctionPlus = wavefunction(rPlus);
            dwavefunction(i,j)=(wavefunctionPlus-wavefunctionMinus)/(2*wavefunctionCurrent*hGrad);
            rPlus(i,j)=r(i,j);
            rMinus(i,j)=r(i,j);
        }
    }
    return dwavefunction;

}


/************************************************************
Name:               loadConfigurations
Description:        loads different variabels
*/
void Wavefunction::loadConfiguration(Config *cfg){
    useAnalyticLaplace=cfg->lookup("AppSettings.useAnalyticLaplace");
    useAnalyticGradient=cfg->lookup("AppSettings.useAnalyticGradient");

}





