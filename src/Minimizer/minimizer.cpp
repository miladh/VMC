#include "minimizer.h"



Minimizer::Minimizer(const int &myRank, const int &nProcess)
{
    this->myRank=myRank;
    this->nProcess=nProcess;
}

/*****************************************************************************
Name:               runMinimizer
Description:        Starts MC-calculation for different variational parameters
*/
void Minimizer::runMinimizaer(){

    alpha=minAlpha;
    beta=minBeta;
    stepAlpha=(maxAlpha-minAlpha)/nVarAlpha;
    stepBeta=(maxBeta-minBeta)/nVarBeta;

    myfile.open ("../vmc/results/results");
    myfile << "Alpha   " << "Beta    " << "Energy       " << "Variance    "<< "Sigma       "<< "Acceptance    "<< endl;

    for(int i=0; i<nVarAlpha;i++){
        vmcapp->alpha=alpha;
        beta=minBeta;

        for(int j=0; j<nVarBeta;j++){

            vmcapp->beta=beta;
            vmcapp->runVMCApp(nCycles,idum);

            Energy=vmcapp->getEnergy();
            EnergySquared=vmcapp->getEnergySquared();
            Variance=vmcapp->getVariance();
            Sigma=vmcapp->getSigma();
            Acceptance= vmcapp->getAcceptanceRate();

            if (myRank == 0) {
                writeToFile();
            }
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
void Minimizer::writeToFile(){
    myfile <<alpha <<"     "<<  beta <<"     "<<Energy
          <<"     "<<Variance<<"     "<<Sigma
         <<"     "<<Acceptance<< endl;

    cout << alpha << ", " << beta << " Energy = " << Energy
         << ", Variance = " << Variance
         << ", Accepted = " << Acceptance
         << "\n";



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
    vmcapp= new VMCApp(cfg,myRank,nProcess);
}


