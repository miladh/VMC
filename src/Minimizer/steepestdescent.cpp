#include <src/Minimizer/steepestdescent.h>

SteepestDescent::SteepestDescent(const int &myRank, const int &nProcess):
    Minimizer(myRank,nProcess),
    nProcess(nProcess),
    myRank(myRank),
    variationalDerivate(zeros<vec>(2))
{
}

/*****************************************************************************
Name:               runMinimizer
Description:        Starts MC-calculation for different variational parameters
*/
void SteepestDescent::runMinimizer(){

    double step = 0.1;
    vec old = zeros<vec>(2);
    for(int i=0; i <100; i++){
        vmcapp->alpha = alpha;
        vmcapp->beta  = beta;
        vmcapp->runVMCApp(nCycles,idum);
        vmcapp->messagePassing();

        variationalDerivate = vmcapp->getVariationalDerivate();
        energy = vmcapp->getEnergy();
        energySquared = vmcapp->getEnergySquared();
        variance=vmcapp->getVariance();
        sigma=vmcapp->getSigma();
        acceptance= vmcapp->getAcceptanceRate();


        if(old(0)/variationalDerivate(0) >= 0){
            step *=1.25;
            alpha -= step*signFunc(variationalDerivate(0));
        }else{
            step*=0.5;
            alpha -= step*signFunc(variationalDerivate(0));
        }


        if(old(1)/variationalDerivate(1) >= 0){
            step*=1.25;
            beta -= step*signFunc(variationalDerivate(0));
        }else{
            step*=0.5;
            beta -= step*signFunc(variationalDerivate(0));
        }

        old = variationalDerivate;


//        if(myRank==0){
//            cout << alpha << ", " << beta << " Energy = " << energy
//                 << ", Variance = " << variance
//                 << ", Accepted = " << acceptance
//                 << "\n";
//        }

    }

}


int SteepestDescent::signFunc(double varDer){

    if(varDer < 0){
        return -1;
    }else{
        return 1;
    }
}



/************************************************************
Name:               loadConfigurations
Description:        loads different variabels
*/
void SteepestDescent::loadConfiguration(Config *cfg){
    alpha=cfg->lookup("setup.MinimizerSettings.SDMinSettings.InitAlpha");
    beta=cfg->lookup("setup.MinimizerSettings.SDMinSettings.InitBeta");
    nCycles=cfg->lookup("AppSettings.cycles");
    idum = cfg->lookup("AppSettings.idum");
    vmcapp= new VMCApp(myRank,nProcess);
    vmcapp->loadConfiguration(cfg);
}
