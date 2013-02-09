#include "minimizer.h"

Minimizer::Minimizer():
    minAlpha(1.6),
    maxAlpha(2.2),
    minBeta(0.2),
    maxBeta(0.4),
    nVarAlpha(2),
    nVarBeta(2),

    Energy(zeros<mat>(nVarAlpha,nVarBeta)),
    EnergySquared(zeros<mat>(nVarAlpha,nVarBeta))
{
}


void Minimizer::runMinimizaer(){

    alpha=minAlpha;
    beta=minBeta;
    stepAlpha=(maxAlpha-minAlpha)/nVarAlpha;
    stepBeta=(maxBeta-minBeta)/nVarBeta;

    for(int i=0; i<nVarAlpha;i++){
        vmcapp.alpha=alpha;
        beta=minBeta;

        for(int j=0; j<nVarBeta;j++){

            vmcapp.beta=beta;
            vmcapp.runVMCApp();

            Energy(i,j)=vmcapp.energy;
            EnergySquared(i,j)=vmcapp.energySquared;

            cout << "alpha: " <<alpha << " beta: "<< beta<< " Energy: "<< Energy(i,j) <<endl;

            beta+=stepBeta;

        }
        alpha+=stepAlpha;
    }

}
