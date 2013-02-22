#include "wavefunction.h"

Wavefunction::Wavefunction()
{
    orbitals = new Hydrogenic;
    hGrad= 1e-5;
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


    hVec=h*ones<rowvec>(r.n_cols);
    rPlus = zeros<mat>( r.n_cols);
    rMinus = zeros<mat>(r.n_cols);

    rPlus = rMinus = r;

    wavefunctionMinus = 0;
    wavefunctionPlus = 0;

    wavefunctionCurrent = wavefunction(r);


    for(uint i = 0; i < r.n_rows; i++) {
        rPlus.row(i) += hVec;
        rMinus.row(i) -= hVec;
        wavefunctionMinus = wavefunction(rMinus);
        wavefunctionPlus = wavefunction(rPlus);
        dwavefunction.row(i)= ones<rowvec>(r.n_cols)*(wavefunctionPlus-wavefunctionMinus)/(wavefunctionCurrent*hGrad);
        rPlus.row(i) = r.row(i);
        rMinus.row(i)= r.row(i);
    }

    return dwavefunction;

}


/************************************************************
Name:               loadConfigurations
Description:        loads different variabels
*/
void Wavefunction::loadConfiguration(Config *cfg){
    charge=cfg->lookup("PotentialSettings.charge");
    useAnalyticLaplace=cfg->lookup("AppSettings.useAnalyticLaplace");
    useAnalyticGradient=cfg->lookup("AppSettings.useAnalyticGradient");

}





