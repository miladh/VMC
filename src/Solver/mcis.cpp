#include "mcis.h"
#include "src/includes/lib.h"
#include "src/Hamiltonian/hamiltonian.h"
#include "src/Potential/coulomb_potential.h"
#include "src/Kinetic/kinetic.h"
#include <armadillo>
#include <iostream>
#include <math.h>
#include <iomanip>

MCIS::MCIS(const uint &nParticles, const uint &nDimensions,Hamiltonian *hamiltonian, Wavefunction* TrialWavefunction):
    Solver(nParticles,nDimensions,hamiltonian,TrialWavefunction),
    qForceOld (zeros<mat>(nParticles, nDimensions)),
    qForceNew (zeros<mat>(nParticles, nDimensions)),
    timeStep(0.005),
    D(0.5)
{
}



/************************************************************
Name:               solve
Description:        starts a MC-sample
*/
void MCIS::solve(int nCycles, long idum)
{
    this->idum=idum;
    this->nCycles=nCycles;
    MetropolisAlgoIS();
}



/************************************************************
Name:               MetropolisAlgoIS
Description:        MetropolisAlgo important sampling
*/
void MCIS::MetropolisAlgoIS(){
    acceptedSteps=0;
    energySum = 0;
    energySquaredSum = 0;

    rOld = randn(nParticles,nDimensions)*sqrt(timeStep); //std=sqrt(2*D*dt), D=0.5
    rNew = rOld;
    TrialWavefunction->initializewavefunction(rOld);
    qForceOld = getQuantumForce(rOld);

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles+thermalization; cycle++) {
        // New position to test
        for(int i = 0; i < nParticles; i++) {

            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + randn()*sqrt(timeStep)+qForceOld(i,j)*timeStep*D;
            }

            TrialWavefunction->activeParticle(rNew,i);
            TrialWavefunction->updateWavefunction();

            qForceNew=getQuantumForce(rNew);

            //compute green's function ratio
            GreensFunction=0;
            for (int k = 0; k < nParticles; k++){
                for(int l=0; l<nDimensions; l++){
                    GreensFunction+=(qForceOld(k,l)+qForceNew(k,l))*
                            (D*timeStep*0.5*(qForceOld(k,l)-qForceNew(k,l))-rNew(k,l)+rOld(k,l));
                }
            }

            GreensFunction= exp(0.5*GreensFunction);


            R =TrialWavefunction->getRatio();
            R*=R*GreensFunction;

            if(ran2(&idum) <= R) {
                TrialWavefunction->acceptMove();
                rOld.row(i) = rNew.row(i);
                qForceOld=qForceNew;

                if(cycle > thermalization){
                    acceptedSteps++;
                }

            } else {
                TrialWavefunction->rejectMove();
                rNew.row(i) = rOld.row(i);
                qForceNew=qForceOld;
            }

        }
        // update energies
        if(cycle > thermalization){
            deltaE =hamiltonian->getEnergy(rNew);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
    }

    energy = energySum/(nCycles-1);
    energySquared = energySquaredSum/(nCycles-1);
    acceptedSteps= acceptedSteps/nParticles;

}


/************************************************************
Name:               QuantumForce
Description:
*/
mat MCIS::getQuantumForce(const mat &r){

    return 2*TrialWavefunction->gradient(r);
}

