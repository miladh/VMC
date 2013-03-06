#include "mcis.h"
#include "src/includes/lib.h"
#include "src/Hamiltonian/hamiltonian.h"
#include "src/Potential/coulomb_potential.h"
#include "src/Kinetic/kinetic.h"
#include <armadillo>
#include <iostream>
#include <math.h>
#include <iomanip>

MCIS::MCIS(Hamiltonian *hamiltonian, Wavefunction* TrialWavefunction)
{
    this->hamiltonian=hamiltonian;
    this->TrialWavefunction=TrialWavefunction;
    timeStep=0.02;
    D=0.5;
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
    qForceOld = zeros<mat>(nParticles, nDimensions);
    qForceNew = zeros<mat>(nParticles, nDimensions);
//    qForce = zeros<mat>(nParticles, nDimensions);

    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    wavefunctionOld = 0;
    wavefunctionNew = 0;

    energySum = 0;
    energySquaredSum = 0;


    rOld = randn(nParticles,nDimensions)*sqrt(timeStep); //std=sqrt(2*D*dt), D=0.5
    rNew = rOld;


    wavefunctionOld=TrialWavefunction->wavefunction(rOld);
    qForceOld = getQuantumForce(rOld);


    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles+thermalization; cycle++) {
        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + randn()*sqrt(timeStep)+qForceOld(i,j)*timeStep*D;
            }

            // Recalculate the value of the wave function
            wavefunctionNew = TrialWavefunction->wavefunction(rNew);
            qForceNew=getQuantumForce(rNew);

            //compute green's function ratio
            GreensFunction=0;
            for(int l=0; l<nDimensions; l++){
                GreensFunction+=(qForceOld(i,l)+qForceNew(i,l))*
                        (D*timeStep*0.5*(qForceOld(i,l)-qForceNew(i,l))-rNew(i,l)+rOld(i,l));
            }
            GreensFunction= exp(0.5*GreensFunction);


            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (GreensFunction*wavefunctionNew*wavefunctionNew) / (wavefunctionOld*wavefunctionOld)) {
                rOld.row(i) = rNew.row(i);
                qForceOld.row(i)=qForceNew.row(i);
                wavefunctionOld = wavefunctionNew;

                if(cycle > thermalization){
                    acceptedSteps++;
                }

            } else {
                rNew.row(i) = rOld.row(i);
                qForceNew.row(i)=qForceOld.row(i);

            }

            // update energies
            if(cycle > thermalization){
                deltaE =hamiltonian->getEnergy(rNew);
                energySum += deltaE;
                energySquaredSum += deltaE*deltaE;
            }
        }
    }

    energy = energySum/(nCycles * nParticles);
    energySquared = energySquaredSum/(nCycles * nParticles);
    acceptedSteps= acceptedSteps/nParticles;

}


/************************************************************
Name:               QuantumForce
Description:
*/
mat MCIS::getQuantumForce(const mat &r){

    return 2*TrialWavefunction->gradient(r);
}

