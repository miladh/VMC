#include "mcbf.h"
#include "src/includes/lib.h"
#include "src/Wavefunction/jastrowwavefunction.h"
#include "src/Wavefunction/basicwavefunction.h"
#include "src/Hamiltonian/hamiltonian.h"
#include "src/Potential/coulomb_potential.h"
#include "src/Kinetic/kinetic.h"
#include <armadillo>
#include <iostream>
#include <math.h>
#include <iomanip>

MCBF::MCBF(Hamiltonian *hamiltonian, Wavefunction* TrialWaveFunction)
{
    this->hamiltonian=hamiltonian;
    this->TrialWaveFunction=TrialWaveFunction;
}


/************************************************************
Name:               solve
Description:        starts a MC-sample
*/

void MCBF::solve(int nCycles, long idum)

{
    this->idum=idum;
    stepLength=optimalStepLength();
    MetropolisAlgoBF(nCycles,stepLength);


}


/************************************************************
Name:               MetropolisAlgoBF
Description:
*/
void MCBF::MetropolisAlgoBF(int nCycles, double stepLength){

    acceptedSteps=0;
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    waveFunctionOld = 0;
    waveFunctionNew = 0;

    energySum = 0;
    energySquaredSum = 0;


    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles+thermalization; cycle++) {

        // Store the current value of the wave function
        waveFunctionOld = TrialWaveFunction->wavefunction(rOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }
           // rNew.row(i) =rOld.row(i)+sqrt(stepLength)*randn<rowvec>(nDimensions);

            // Recalculate the value of the wave function
            waveFunctionNew = TrialWaveFunction->wavefunction(rNew);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                    rOld.row(i) = rNew.row(i);
                    waveFunctionOld = waveFunctionNew;

                if(cycle > thermalization){
                    acceptedSteps++;
                }

            } else {
                    rNew.row(i) = rOld.row(i);

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
Name:               optimalStepLength
Description:        Finds the optimal steplength using intersection method
*/

double MCBF::optimalStepLength() {

    while ((maxStepLength - minStepLength) > tolerance) {
        MetropolisAlgoBF(nPreCycles,minStepLength);
        stepMin=acceptedSteps/nPreCycles-0.5;

        MetropolisAlgoBF(nPreCycles,(minStepLength + maxStepLength)/2);
        stepMinMax=acceptedSteps/nPreCycles-0.5;


        if (stepMin*stepMinMax < 0)
            maxStepLength = (minStepLength + maxStepLength) / 2;
        else
            minStepLength = (minStepLength + maxStepLength) / 2;
    }


    return (minStepLength + maxStepLength) / 2;
}




