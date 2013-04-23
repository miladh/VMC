#include <src/Minimizer/bfminimizer.h>


BFMinimizer::BFMinimizer(const int& myRank, const int& nProcess):
    Minimizer(myRank,nProcess)
{
    this->myRank=myRank;
    this->nProcess=nProcess;
}

/*****************************************************************************
Name:               runMinimizer
Description:        Starts MC-calculation for different variational parameters
*/
void BFMinimizer::runMinimizer(){

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
            vmcapp->messagePassing();

            Energy =vmcapp->getEnergy();
            EnergySquared = vmcapp->getEnergySquared();
            Variance = vmcapp->getVariance();
            Sigma = vmcapp->getSigma();
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
void BFMinimizer::writeToFile(){
    myfile <<alpha <<"     "<<  beta <<"     "<<Energy
          <<"     "<<Variance<<"     "<<Sigma
         <<"     "<<Acceptance<< endl;

//    cout << alpha << ", " << beta << " Energy = " << Energy
//         << ", Variance = " << Variance
//         << ", Accepted = " << Acceptance
//         << "\n";



}


/************************************************************
Name:               loadConfigurations
Description:        loads different variabels
*/
void BFMinimizer::loadConfiguration(Config *cfg){
    minAlpha=cfg->lookup("setup.MinimizerSettings.BFMinSettings.minalpha");
    maxAlpha=cfg->lookup("setup.MinimizerSettings.BFMinSettings.maxalpha");
    minBeta=cfg->lookup("setup.MinimizerSettings.BFMinSettings.minbeta");
    maxBeta=cfg->lookup("setup.MinimizerSettings.BFMinSettings.maxbeta");
    nVarAlpha=cfg->lookup("setup.MinimizerSettings.BFMinSettings.nVarAlpha");
    nVarBeta=cfg->lookup("setup.MinimizerSettings.BFMinSettings.nVarBeta");
    nCycles=cfg->lookup("AppSettings.cycles");
    idum=cfg->lookup("AppSettings.idum");

    vmcapp= new VMCApp(myRank,nProcess);
    vmcapp->loadConfiguration(cfg);
}


