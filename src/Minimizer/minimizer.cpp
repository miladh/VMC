#include "minimizer.h"
//#include <mpi.h>

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

    alpha=minAlpha;
    beta=minBeta;
    stepAlpha=(maxAlpha-minAlpha)/nVarAlpha;
    stepBeta=(maxBeta-minBeta)/nVarBeta;

    myfile.open ("../vmc/src/results");
   //myfile << "Alpha   " << "Beta    " << "Energy    " << endl;

    for(int i=0; i<nVarAlpha;i++){
        vmcapp->alpha=alpha;
        beta=minBeta;

        for(int j=0; j<nVarBeta;j++){

            vmcapp->beta=beta;
            vmcapp->runVMCApp();

            Energy(i,j)=vmcapp->energy;
            EnergySquared(i,j)=vmcapp->energySquared;

            cout << "alpha: " <<alpha << "  beta:    "<< beta << "  Energy: "<< Energy(i,j) <<endl;

            myfile <<alpha <<"     "<<  beta <<"     "<<Energy(i,j) <<endl;

            beta+=stepBeta;

        }
        alpha+=stepAlpha;
    }
    myfile.close();
}



/************************************************************
Name:               eloadConfigurations
Description:        loads different variabels
*/
void Minimizer::loadConfiguration(Config *cfg){
    minAlpha=cfg->lookup("MinimizerSettings.minalpha");
    maxAlpha=cfg->lookup("MinimizerSettings.maxalpha");
    minBeta=cfg->lookup("MinimizerSettings.minbeta");
    maxBeta=cfg->lookup("MinimizerSettings.maxbeta");
    nVarAlpha=cfg->lookup("MinimizerSettings.nVarAlpha");
    nVarBeta=cfg->lookup("MinimizerSettings.nVarBeta");
    vmcapp= new VMCApp(cfg);
}
