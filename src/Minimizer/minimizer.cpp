#include "minimizer.h"



Minimizer::Minimizer()
{
}

/*****************************************************************************
Name:               runMinimizer
Description:        Starts MC-calculation for different variational parameters
*/
void Minimizer::runMinimizaer(){

    Energy=zeros<mat>(nVarAlpha,nVarBeta);
    EnergySquared=zeros<mat>(nVarAlpha,nVarBeta);
    Variance=zeros<mat>(nVarAlpha,nVarBeta);
    Acceptance=zeros<mat>(nVarAlpha,nVarBeta);
    Sigma=zeros<mat>(nVarAlpha,nVarBeta);

    alpha=minAlpha;
    beta=minBeta;
    stepAlpha=(maxAlpha-minAlpha)/nVarAlpha;
    stepBeta=(maxBeta-minBeta)/nVarBeta;

    myfile.open ("../vmc/results/results");
   // myfile << "Alpha   " << "Beta    " << "Energy       " << "Variance    "<< "Sigma       "<< "Acceptance    "<< endl;

    for(int i=0; i<nVarAlpha;i++){
        vmcapp->alpha=alpha;
        beta=minBeta;

        for(int j=0; j<nVarBeta;j++){

            vmcapp->beta=beta;
            vmcapp->runVMCApp(nCycles,idum);

            Energy(i,j)=vmcapp->energy;
            EnergySquared(i,j)=vmcapp->energySquared;
//            Variance(i,j)=vmcapp->Variance;
//            Sigma(i,j)=vmcapp->Sigma;
//            Acceptance(i,j)= vmcapp->Acceptance;


            myfile <<alpha <<"     "<<  beta <<"     "<<Energy(i,j)<<endl;

            beta+=stepBeta;

        }
        alpha+=stepAlpha;
    }
    myfile.close();
}



/************************************************************
Name:               loadConfigurations
Description:        loads different variabels
*/
void Minimizer::loadConfiguration(Config *cfg){
    minAlpha=cfg->lookup("MinimizerSettings.minalpha");
    maxAlpha=cfg->lookup("MinimizerSettings.maxalpha");
    minBeta=cfg->lookup("MinimizerSettings.minbeta");
    maxBeta=cfg->lookup("MinimizerSettings.maxbeta");
    nVarAlpha=cfg->lookup("MinimizerSettings.nVarAlpha");
    nVarBeta=cfg->lookup("MinimizerSettings.nVarBeta");
    nCycles=cfg->lookup("AppSettings.cycles");
    idum=cfg->lookup("AppSettings.idum");
    vmcapp= new VMCApp(cfg);
}



