#include "vmcapp.h"
#include <src/includes/Defines.h>
#include <src/Solver/mcbf.h>
#include <src/Solver/mcis.h>
#include <src/includes/lib.h>
#include <src/Wavefunction/hlikewavefunction.h>
#include <src/Jastrow/padejastrow.h>
#include <src/Jastrow/nojastrow.h>
#include <src/Potential/coulombPotential.h>
#include <src/Kinetic/kinetic.h>
#include <src/electronInteraction/coulombinteraction.h>
#include <src/electronInteraction/nointeraction.h>



VMCApp::VMCApp(const int &myRank, const int &nProcess):
    nProcess(nProcess),
    myRank(myRank)
{
}

/************************************************************
Name:               runVMCApp
Description:        starts VMC calculations
*/
void VMCApp::runVMCApp(int nCycles, long idum)
{

    idum = idum - myRank - time(NULL);
    srand(-idum);

    nCycles /= nProcess;

    trialWavefunction = setWavefunction();
    trialWavefunction->loadConfiguration(cfg);


    potential = new CoulombPotential;
    potential->loadConfiguration(cfg);

    kinetic= new Kinetic;
    kinetic->wf=trialWavefunction;

    electonInteraction = setInteraction();

    hamiltonian =new Hamiltonian;
    hamiltonian->potential=potential;
    hamiltonian->kinetic=kinetic;
    hamiltonian->electronInteraction = electonInteraction;

    observables = new Observables(hamiltonian,trialWavefunction);
    observables->loadConfiguration(cfg);

    solver = setSolverMethod();
    solver->loadConfiguration(cfg);
    solver->initializeSolver();
    solver->solve(nCycles,idum);


}
/************************************************************
Name:
Description:
*/
void VMCApp::messagePassing()
{

    tmp = observables->getEnergy();
    MPI_Allreduce(&tmp, &totEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    totEnergy/= nProcess;

    tmp = observables->getEnergySquared();
    MPI_Allreduce(&tmp, &totEnergySquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    totEnergySquared /= nProcess;

    tmp = solver->acceptedSteps;
    MPI_Allreduce(&tmp, &Acceptance, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Acceptance /= nProcess;

    if(minimize){
        tmpVec = observables->getVariationalDerivateRatio();
        MPI_Allreduce(&tmpVec[0], &totVariationalDerivate[0], tmpVec.n_rows, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        totVariationalDerivate /= nProcess;

        tmpVec = observables->getEnergyVariationalDerivate();
        MPI_Allreduce(&tmpVec[0], &totEnergyVarDerivate[0], tmpVec.n_rows, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        totEnergyVarDerivate /= nProcess;
    }

    if(doBlocking){
        observables->writeEnergyVectorToFile(myRank);
    }
}


/************************************************************
Name:               setSolverMethod
Description:
*/
ElectronInteraction *VMCApp::setInteraction(){

    ElectronInteraction* interaction;
    InteractionType = cfg->lookup("AppSettings.interactionType");

    switch (InteractionType) {
    case NOINTERACTION:
        interaction = new NoInteraction;
        break;

    case COULOMBINTERACTION:
        interaction =new CoulombInteraction;
        break;
    }
    return interaction;
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
        solver = new MCBF(hamiltonian,trialWavefunction,observables);
        break;
    case IS:
        solver =new MCIS(hamiltonian,trialWavefunction,observables);
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
        wf = new HLikeWavefunction(nParticles);
        wf->jas=new PadeJastrow(nParticles);
        wf->jas->alpha=alpha;
        wf->jas->beta=beta;
        wf->jas->setaValues(nParticles);
        wf->slater->orbitals->k=alpha;
        break;

    case Basic:
        wf = new HLikeWavefunction(nParticles);
        wf->jas = new NoJastrow(nParticles);
        wf->slater->orbitals->k=alpha;
        break;

    case  Hydrogenic:
        wf = new HLikeWavefunction(nParticles);
        wf->jas=new NoJastrow(nParticles);
        wf->slater->orbitals->k=charge;
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

/************************************************************
Name:               getAcceptanceRate
Description:
*/
vec VMCApp::getVariationalDerivate(){
    return 2*totEnergyVarDerivate - 2*totVariationalDerivate*totEnergy;
}


/************************************************************
Name:               loadConfiguration
Description:        loads different variables
*/
void VMCApp::loadConfiguration(Config *cfg){
    this->cfg=cfg;
    nParticles = cfg->lookup("SolverSettings.N");
    nDimensions = cfg->lookup("SolverSettings.dim");
    charge = cfg->lookup("PotentialSettings.charge");
    minimize=cfg->lookup("MinimizerSettings.minimize");
    doBlocking= cfg->lookup("BlockingSettings.doBlocking");
    totVariationalDerivate=zeros<vec>(2);
    totEnergyVarDerivate=zeros<vec>(2);
}
