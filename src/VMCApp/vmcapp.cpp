#include "vmcapp.h"
#include "src/includes/Defines.h"
#include "src/Solver/mcbf.h"
#include "src/Solver/mcis.h"
#include "src/includes/lib.h"
#include "src/Wavefunction/jastrowwavefunction.h"
#include "src/Wavefunction/basicwavefunction.h"
#include "src/Wavefunction/hydrogenicwavefunction.h"
#include "src/Potential/coulomb_potential.h"
#include "src/Kinetic/kinetic.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
#pragma GCC diagnostic warning "-Wunused-parameter"



VMCApp::VMCApp(Config *cfg, const int &myRank, const int &nProcess)
{
    this->cfg=cfg;
    nParticles=cfg->lookup("SolverSettings.N");
    charge=cfg->lookup("PotentialSettings.charge");
    this->myRank=myRank;
    this->nProcess=nProcess;
}

/************************************************************
Name:               runVMCApp
Description:        starts VMC calculations
*/
void VMCApp::runVMCApp(int nCycles, long idum)
{

    idum = idum- myRank - time(NULL);
    nCycles /= nProcess;

    TrialWavefunction = setWavefunction();
    TrialWavefunction->loadConfiguration(cfg);


    potential = new CoulombPotential;
    potential->loadConfiguration(cfg);

    kinetic= new Kinetic;
    kinetic->wf=TrialWavefunction;

    hamiltonian =new Hamiltonian;
    hamiltonian->potential=potential;
    hamiltonian->kinetic=kinetic;

    solver = setSolverMethod();
    solver->loadConfiguration(cfg);
    solver->solve(nCycles,idum);


    tmp = solver->energy;
    MPI_Allreduce(&tmp, &totEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    totEnergy/= nProcess;

    tmp = solver->energySquared;
    MPI_Allreduce(&tmp, &totEnergySquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    totEnergySquared /= nProcess;


    tmp = solver->acceptedSteps;
    MPI_Allreduce(&tmp, &Acceptance, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Acceptance /= (nProcess*nCycles);

}


/************************************************************
Name:               setSolverMethod
Description:
*/
Solver* VMCApp::setSolverMethod(){

    Solver* solver;
    solverType= cfg->lookup("AppSettings.solverType");

    switch (solverType) {
    case BF:
        solver = new MCBF(hamiltonian,TrialWavefunction);
        break;
    case IS:
        solver =new MCIS(hamiltonian,TrialWavefunction);
        break;
    }
    return solver;
}


/************************************************************
Name:               setWaveFunction
Description:
*/
Wavefunction* VMCApp::setWavefunction(){

    Wavefunction* wf;
    WavefunctionType= cfg->lookup("AppSettings.wavefunction");

    switch (WavefunctionType) {
    case Jastrow:
        wf = new JastrowWavefunction;
        wf->jas.alpha=alpha;
        wf->jas.beta=beta;
        wf->jas.setaValues(nParticles);
        wf->orbitals->k=alpha;
        break;

    case Basic:
        wf =new BasicWavefunction;
        wf->jas.alpha=alpha;
        wf->orbitals->k=alpha;
        break;

    case  Hydrogenic:
        wf = new HydrogenicWavefunction(charge);
        wf->orbitals->k=charge;
        break;
    }

    return wf;
}


/************************************************************
Name:               setWaveFunction
Description:
*/
double VMCApp::getEnergy(){

    return totEnergy;
}


/************************************************************
Name:               getEnergySquared
Description:
*/
double VMCApp::getEnergySquared(){

    return totEnergySquared;
}

/************************************************************
Name:               getVariance
Description:
*/
double VMCApp::getVariance(){

    Variance = totEnergySquared - totEnergy * totEnergy;
    return Variance;
}

/************************************************************
Name:               getSigma
Description:
*/
double VMCApp::getSigma(){

    Sigma = sqrt(Variance);
    return Sigma;
}

/************************************************************
Name:               getAcceptanceRate
Description:
*/
double VMCApp::getAcceptanceRate(){

    return Acceptance;
}
