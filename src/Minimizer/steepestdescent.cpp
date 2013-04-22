#include <src/Minimizer/steepestdescent.h>

SteepestDescent::SteepestDescent(const int &myRank, const int &nProcess):
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

    vec old = zeros<vec>(2);
    for(int i=0; i <100; i++){
        vmcapp->alpha = alpha;
        vmcapp->beta  = beta;
        vmcapp->runVMCApp(nCycles,idum);

        variationalDerivate = vmcapp->getVariationalDerivate();
        energy = vmcapp->getEnergy();
        energySquared = vmcapp->getEnergySquared();
        variance=vmcapp->getVariance();
        sigma=vmcapp->getSigma();
        acceptance= vmcapp->getAcceptanceRate();


        if(old(0)/variationalDerivate(0) >= 0){
            alpha -= variationalDerivate(0)*1.1;
        }else{
            alpha -= variationalDerivate(0)*0.5;
        }


        if(old(1)/variationalDerivate(1) >= 0){
            beta -= variationalDerivate(1)*1.1;
        }else{
            beta -= variationalDerivate(1)*0.5;
        }

        old = variationalDerivate(0);


        if(myRank==0){
            cout << alpha << ", " << beta << " Energy = " << energy
                 << ", Variance = " << variance
                 << ", Accepted = " << acceptance
                 << "\n";
        }

    }

}


/************************************************************
Name:               loadConfigurations
Description:        loads different variabels
*/
void SteepestDescent::loadConfiguration(Config *cfg){
    alpha=cfg->lookup("MinimizerSettings.minalpha");
    beta=cfg->lookup("MinimizerSettings.minbeta");
    nCycles=cfg->lookup("AppSettings.cycles");
    idum=cfg->lookup("AppSettings.idum");
    vmcapp= new VMCApp(cfg,myRank,nProcess);
}
