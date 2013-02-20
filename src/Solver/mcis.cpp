#include "mcis.h"
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

MCIS::MCIS(Hamiltonian *hamiltonian, Wavefunction* TrialWaveFunction)
{
    this->hamiltonian=hamiltonian;
    this->TrialWaveFunction=TrialWaveFunction;
}



/************************************************************
Name:               solve
Description:        starts a MC-sample
*/
void MCIS::solve(int nCycles, long idum)

{
    this->idum=idum;
    timeStep=0.02;
    D=0.5;
    MetropolisAlgoIS(nCycles,idum);

}

/************************************************************
Name:               MetropolisAlgoIS
Description:        MetropolisAlgo important sampling
*/
void MCIS::MetropolisAlgoIS(int nCycles,long idum){
    acceptedSteps=0;

    qForceOld = zeros<mat>(nParticles, nDimensions);
    qForceNew = zeros<mat>(nParticles, nDimensions);
    qForce = zeros<mat>(nParticles, nDimensions);


    waveFunctionOld = 0;
    waveFunctionNew = 0;

    energySum = 0;
    energySquaredSum = 0;




    rOld = randn(nParticles,nDimensions)*sqrt(timeStep);
    rNew = rOld;


    waveFunctionOld=TrialWaveFunction->waveFunction(nParticles,rOld);
    qForceOld=QuantumForce(rOld);




    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles+thermalization; cycle++) {

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + randn()*sqrt(timeStep)+qForceOld(i,j)*timeStep*D;
            }

            // Recalculate the value of the wave function
            waveFunctionNew = TrialWaveFunction->waveFunction(nParticles,rNew);
            qForceNew=QuantumForce(rNew);


            GreensFunction=0;
            for(int l=0; l<nDimensions; l++){
                GreensFunction+=(qForceOld(i,l)+qForceNew(i,l))*(D*timeStep*0.5*(qForceOld(i,l)-qForceNew(i,l))-rNew(i,l)+rOld(i,l));
            }
            GreensFunction= exp(0.5*GreensFunction);


            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (GreensFunction*waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                rOld.row(i) = rNew.row(i);
                qForceOld.row(i)=qForceNew.row(i);
                waveFunctionOld = waveFunctionNew;

                if(cycle > thermalization){
                    acceptedSteps++;
                }

            } else {
                rNew.row(i) = rOld.row(i);
                qForceNew.row(i)=qForceOld.row(i);

            }

            // update energies
            if(cycle > thermalization){
                deltaE =hamiltonian->getEnergy(nParticles,rNew);
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
mat MCIS::QuantumForce(const mat &r){

    rowvec hVec=h*ones<rowvec>(r.n_cols);
    mat rPlus = zeros<mat>(nParticles, r.n_cols);
    mat rMinus = zeros<mat>(nParticles, r.n_cols);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = TrialWaveFunction->waveFunction(nParticles,r);

    for(int i = 0; i < nParticles; i++) {
        rPlus.row(i) += hVec;
        rMinus.row(i) -= hVec;
        waveFunctionMinus = TrialWaveFunction->waveFunction(nParticles,rMinus);
        waveFunctionPlus = TrialWaveFunction->waveFunction(nParticles,rPlus);
        qForce.row(i)= ones<rowvec>(r.n_cols)*(waveFunctionPlus-waveFunctionMinus)/(waveFunctionCurrent*h);
        rPlus.row(i) = r.row(i);
        rMinus.row(i)= r.row(i);
    }

    return qForce;

}


/************************************************************
Name:               loadConfiguration
Description:        loads different variables
*/
void MCIS::loadConfiguration(Config *cfg){
    nDimensions=cfg->lookup("SolverSettings.dim");
    nParticles=cfg->lookup("SolverSettings.N");
    thermalization=cfg->lookup("AppSettings.thermalization");
    h = cfg->lookup("NumericalKineticSettings.h");

}



//            for(int k=0; k <nParticles; k++){
//                if(k!=i){
//                    rNew.row(k) = rOld.row(k);
//                }
//            }
