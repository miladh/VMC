#include "mcis.h"



MCIS::MCIS(Hamiltonian *hamiltonian, Wavefunction* TrialWavefunction, Observables* observables):
    Solver(hamiltonian,TrialWavefunction,observables),
    timeStep(0.05),
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
    vector <mat> positionsMat;
    observables->initializeObservables(nCycles);
    qForceOld = zeros<mat>(nParticles, nDimensions);
    qForceNew = zeros<mat>(nParticles, nDimensions);
    acceptedSteps=0;


    rOld = randn(nParticles,nDimensions)*sqrt(timeStep); //std=sqrt(2*D*dt), D=0.5
    rNew = rOld;
    trialWavefunction->initializewavefunction(rOld);
    qForceOld = getQuantumForce(rOld);

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles+thermalization; cycle++) {
        // New position to test
        for(int i = 0; i < nParticles; i++) {

            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + randn()*sqrt(timeStep)+qForceOld(i,j)*timeStep*D;
            }

            trialWavefunction->activeParticle(rNew,i);
            trialWavefunction->updateWavefunction();

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


            R =trialWavefunction->getRatio();
            R*=R*GreensFunction;

            if(ran2(&idum) <= R) {
                trialWavefunction->acceptMove();
                rOld.row(i) = rNew.row(i);
                qForceOld=qForceNew;

                if(cycle > thermalization){
                    acceptedSteps++;

//                    if(cycle%25==0){
//                        positionsMat.push_back(rNew);
//                    }

                }

            } else {
                trialWavefunction->rejectMove();
                rNew.row(i) = rOld.row(i);
                qForceNew=qForceOld;
            }

        }
        // update energies
        if(cycle > thermalization){
            observables->currentConfiguration(rNew);
            observables->calculateObservables();
        }
    }

    //    if(myRank==0){
    //        ofstream myfile;
    //        myfile.open("../vmc/results/onebodyDensity/OBD.mat");
    //        for(uint i=0; i<positionsMat.size(); i++){
    //            myfile << positionsMat[i] <<endl;
    //        }
    //    }
        acceptedSteps /= ((nCycles-1)*nParticles);
}


/************************************************************
Name:               QuantumForce
Description:
*/
mat MCIS::getQuantumForce(const mat &r){

    return 2*trialWavefunction->gradient(r);
}

