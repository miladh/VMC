#include "vmcsolver.h"
#include "includes/lib.h"
#include "src/Wavefunction/jastrowwavefunction.h"
#include "src/Wavefunction/basicwavefunction.h"
#include "src/Potential/coulomb_potential.h"
#include "src/Kinetic/numericalkinetic.h"
#include "src/Kinetic/closedformkinetic.h"
#include <armadillo>
#include <iostream>


#include <typeinfo>

using namespace arma;
using namespace std;

VMCSolver::VMCSolver() :
    nDimensions(3),
    charge(2),
    stepLength(1.0),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(-1),
    alpha(1.85),
    beta(0.25),
    nCycles(1000000)
{
    TrialWaveFunction = new JastrowWaveFunction();
//    TrialWaveFunction = new BasicWaveFunction();
    TrialWaveFunction->alpha=alpha;
    TrialWaveFunction->beta=beta;

    PotentialEnergy = new CoulombPotential();

    KineticEnergy = new ClosedFormKinetic();
    KineticEnergy->wf = TrialWaveFunction;
    KineticEnergy->alpha=alpha;
    KineticEnergy->beta=beta;

}

void VMCSolver::runMonteCarloIntegration()
{
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

    double energySum = 0;
    double energySquaredSum = 0;

    double deltaE;

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
        waveFunctionOld = TrialWaveFunction->waveFunction(nDimensions,nParticles,rOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }

            // Recalculate the value of the wave function
            waveFunctionNew = TrialWaveFunction->waveFunction(nDimensions,nParticles,rNew);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    waveFunctionOld = waveFunctionNew;
                }
            } else {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
            }
            // update energies
            deltaE = localEnergy(rNew);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
}


double VMCSolver::localEnergy(const mat &r)
{
    return KineticEnergy->evaluate(nDimensions,nParticles,r);// + PotentialEnergy->evaluate(nDimensions,nParticles,r);
}


