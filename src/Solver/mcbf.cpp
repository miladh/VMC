#include "mcbf.h"
#include "src/includes/lib.h"
#include "src/Wavefunction/jastrowwavefunction.h"
#include "src/Wavefunction/basicwavefunction.h"
#include "src/Hamiltonian/hamiltonian.h"
#include "src/Potential/coulomb_potential.h"
#include "src/Kinetic/numericalkinetic.h"
#include "src/Kinetic/closedformkinetic.h"
#include <armadillo>
#include <iostream>
#include <math.h>


MCBF::MCBF(Hamiltonian *hamiltonian, Wavefunction* TrialWaveFunction):
    nDimensions(3),
    nParticles(2),
    nCycles(1000000),
    idum(-1)
{
    this->hamiltonian=hamiltonian;
    this->TrialWaveFunction=TrialWaveFunction;

}

void MCBF::solve()

{
    stepLength=optimalStepLength();
    //cout << "step"<<stepLength<<endl;
    MetropolisAlgo(nCycles,stepLength);

}

void MCBF::MetropolisAlgo(int nCycles, double stepLength){
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
    for(int cycle = 0; cycle < nCycles; cycle++) {

        // Store the current value of the wave function
        waveFunctionOld = TrialWaveFunction->waveFunction(nParticles,rOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
                //rNew(i,j) = rOld(i,j) + sqrt(stepLength)*as_scalar(randn(1));
            }
            //rNew.row(i) += sqrt(stepLength)*randn<rowvec>(nDimensions);

            // Recalculate the value of the wave function
            waveFunctionNew = TrialWaveFunction->waveFunction(nParticles,rNew);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    waveFunctionOld = waveFunctionNew;
                }

            acceptedSteps++;

            } else {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
            }
            // update energies
            deltaE =hamiltonian->getEnergy(nParticles,rNew);

            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
    }
    energy = energySum/(nCycles * nParticles);
    energySquared = energySquaredSum/(nCycles * nParticles);
    acceptedSteps= acceptedSteps/nParticles;
    cout << "accepted "<< acceptedSteps/(nCycles)<<endl;

}



double MCBF::optimalStepLength() {
    nPreCycles=10000;
    double minStepLength, maxStepLength, tolerance;
    double stepMinMax,stepMin;
    minStepLength = -3;
    maxStepLength = 3;
    tolerance = 0.01;

    while ((maxStepLength - minStepLength) > tolerance) {
        MetropolisAlgo(nPreCycles,minStepLength);
        stepMin=acceptedSteps/nPreCycles-0.5;

        MetropolisAlgo(nPreCycles,(minStepLength + maxStepLength)/2);
        stepMinMax=acceptedSteps/nPreCycles-0.5;


        if (stepMin*stepMinMax < 0)
            maxStepLength = (minStepLength + maxStepLength) / 2;
        else
            minStepLength = (minStepLength + maxStepLength) / 2;
    }


    return (minStepLength + maxStepLength) / 2;
}



//    nPreCycles=1000;
//    double h=0.01;
//    double xOld=0;
//    double xNew=3;
//    double plus,minus;
//    double f,fminus,fplus,ff;


//    while(abs(xNew-xOld)>=0.01){
//        xOld=xNew;
//        minus=xOld-h;
//        plus=xOld+h;

//        MetropolisAlgo(nPreCycles,xOld);
//        f= acceptedSteps/nPreCycles-0.5;

//        MetropolisAlgo(nPreCycles,minus);
//        fminus= acceptedSteps/nPreCycles-0.5;

//        MetropolisAlgo(nPreCycles,plus);
//        fplus= acceptedSteps/nPreCycles-0.5;

//        ff= (fplus-fminus)/(2*h);
//        xNew= xOld -f/ff;
//    }
//    return xNew;

//}


