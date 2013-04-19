#include "mcbf.h"
#include <src/includes/lib.h>
#include <src/Hamiltonian/hamiltonian.h>
#include <src/Potential/coulombPotential.h>
#include <src/Kinetic/kinetic.h>


MCBF::MCBF(const uint &nParticles, const uint &nDimensions, Hamiltonian *hamiltonian, Wavefunction* TrialWavefunction):
    Solver(nParticles,nDimensions,hamiltonian,TrialWavefunction)
{
}


/************************************************************
Name:               solve
Description:        starts a MC-sample
*/

void MCBF::solve(int nCycles, long idum)
{
    this->idum=idum;

#if BLOCKING
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    energyVector=zeros(nCycles-1);
#endif

    stepLength=optimalStepLength();
    MetropolisAlgoBF(nCycles,stepLength);

}


/************************************************************
Name:               MetropolisAlgoBF
Description:
*/
void MCBF::MetropolisAlgoBF(int nCycles, double stepLength){

    acceptedSteps=0;
    energySum = 0;
    energySquaredSum = 0;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    // Store the current value of the wave function
    TrialWavefunction->initializewavefunction(rOld);



    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles+thermalization; cycle++) {

        // New position to test
        for(int i = 0; i < nParticles; i++) {

            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }

            // Recalculate the value of the wave function
            TrialWavefunction->activeParticle(rNew,i);
            TrialWavefunction->updateWavefunction();
            R =TrialWavefunction->getRatio();
            R*=R;

            if(ran2(&idum) <= R) {
                TrialWavefunction->acceptMove();
                rOld.row(i) = rNew.row(i);

                if(cycle > thermalization){
                    acceptedSteps++;
                }

            } else {
                TrialWavefunction->rejectMove();
                rNew.row(i) = rOld.row(i);

            }

        }

        // update energies
        if(cycle > thermalization){
            deltaE =hamiltonian->getEnergy(rNew);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;       
#if BLOCKING
            energyVector(cycle-thermalization-1)=deltaE;
#endif
        }
    }

    energy = energySum/(nCycles-1);
    energySquared = energySquaredSum/(nCycles-1);
    acceptedSteps= acceptedSteps/nParticles;

#if BLOCKING
    ostringstream filename;
    filename << "../vmc/results/blocking/data_" << myRank << ".mat";
    energyVector.save(filename.str());
#endif
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




