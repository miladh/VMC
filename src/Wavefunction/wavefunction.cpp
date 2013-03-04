#include "wavefunction.h"

Wavefunction::Wavefunction()
{
    orbitals = new Hydrogenic;
    hGrad= 1e-8;
    h =0.001;
    h2 =1000000.0;
}


/************************************************************
Name:               laplaceNumerical
Description:
*/
double Wavefunction::laplaceNumerical(const mat &r){

    hVec=h*ones<rowvec>(r.n_cols);
    rPlus = zeros<mat>(r.n_rows, r.n_cols);
    rMinus = zeros<mat>(r.n_rows, r.n_cols);

    rPlus = rMinus = r;

    wavefunctionMinus = 0;
    wavefunctionPlus = 0;

    wavefunctionCurrent = wavefunction(r);

    ddwavefunction = 0;
    for(uint i = 0; i <r.n_rows; i++) {
        rPlus.row(i) += hVec;
        rMinus.row(i) -= hVec;
        wavefunctionMinus = wavefunction(rMinus);
        wavefunctionPlus = wavefunction(rPlus);
        ddwavefunction += (wavefunctionMinus + wavefunctionPlus - 2 * wavefunctionCurrent);
        rPlus.row(i) = r.row(i);
        rMinus.row(i)= r.row(i);
    }

    ddwavefunction =h2 * ddwavefunction / wavefunctionCurrent;

    return ddwavefunction;
}



/************************************************************
Name:               gradientNumerical
Description:
*/
mat Wavefunction::gradientNumerical(const mat &r){

    rPlus = zeros<mat>(r.n_rows,r.n_cols);
    rMinus = zeros<mat>(r.n_rows,r.n_cols);
    dwavefunction=zeros<mat>(r.n_rows,r.n_cols);
    rPlus = rMinus = r;
    wavefunctionMinus = 0.0;
    wavefunctionPlus = 0.0;

    wavefunctionCurrent = wavefunction(r);


    for(uint i =0; i<r.n_rows ; i++){
        for(uint j =0; j<r.n_cols ; j++){
            rPlus(i,j)+=hGrad;
            rMinus(i,j)-=hGrad;
            wavefunctionMinus = wavefunction(rMinus);
            wavefunctionPlus = wavefunction(rPlus);
            dwavefunction(i,j)=(wavefunctionPlus-wavefunctionMinus)/(wavefunctionCurrent*hGrad);
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





