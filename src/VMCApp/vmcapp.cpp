#include "vmcapp.h"
#include "src/Solver/mcbf.h"
#include "src/Solver/mcis.h"
#include "src/includes/lib.h"
#include "src/Wavefunction/jastrowwavefunction.h"
#include "src/Wavefunction/basicwavefunction.h"
#include "src/Wavefunction/hydrogenicwavefunction.h"
#include "src/Potential/coulomb_potential.h"
#include "src/Kinetic/numericalkinetic.h"
#include "src/Kinetic/closedformkinetic.h"
#include <mpi.h>
#include <iomanip>


VMCApp::VMCApp(Config *cfg)
{
    this->cfg=cfg;
}

/************************************************************
Name:               runVMCApp
Description:        starts VMC calculations
*/
void VMCApp::runVMCApp(int nCycles, long idum)
{

    MPI_Comm_size(MPI_COMM_WORLD, &nProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    idum = idum - myRank - time(NULL);

    nCycles /= nProcess;

    TrialWaveFunction = new JastrowWavefunction;
    TrialWaveFunction->alpha=alpha;
    TrialWaveFunction->beta=beta;

    potential = new CoulombPotential(cfg);

    kinetic= new ClosedFormKinetic(cfg);
    kinetic->wf = TrialWaveFunction;
    kinetic->alpha=alpha;
    kinetic->beta=beta;

    hamiltonian =new Hamiltonian();
    hamiltonian->potential=potential;
    hamiltonian->kinetic=kinetic;

    solver = new MCIS(hamiltonian,TrialWaveFunction);
    solver->loadConfiguration(cfg);
    solver->solve(nCycles,idum);


    tmp = solver->energy;
    MPI_Allreduce(&tmp, &totEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    totEnergy/= nProcess;

    tmp = solver->energySquared;
    MPI_Allreduce(&tmp, &totEnergySquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    totEnergySquared /= nProcess;


    energy= totEnergy;
    energySquared =totEnergySquared;
    Variance = totEnergySquared - totEnergy * totEnergy;
    Sigma = sqrt(Variance);
    Acceptance = solver->acceptedSteps/nCycles;


    if (myRank == 0) {
        cout << alpha << ", " << beta << " Energy = " << totEnergy
             << ", Variance = " << Variance
             //<< ", Sigma = " << sqrt(totEnergySquared - totEnergy * totEnergy)
             << ", Accepted = " << solver->acceptedSteps / nCycles
             //<< ", MC cycles = " << nCycles * nProcess
             << "\n";
    }

}

//void VMCApp::writeToFile(ofstream myfile){
//    this->myfile=myfile;

//    if (myRank == 0) {
//        myfile <<alpha <<"     "<<  beta <<"     "<<energy
//               <<"     "<<Variance <<"     "<<Sigma
//               <<"     "<<Acceptance<< endl;
//    }

//}




