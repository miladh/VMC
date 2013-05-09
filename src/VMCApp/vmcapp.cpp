#include "vmcapp.h"
#include <src/includes/Defines.h>
#include <src/Solver/mcbf.h>
#include <src/Solver/mcis.h>
#include <src/includes/lib.h>
#include <src/Jastrow/padejastrow.h>
#include <src/Jastrow/nojastrow.h>
#include <src/Potential/coulombPotential.h>
#include <src/Kinetic/kinetic.h>
#include <src/electronInteraction/coulombinteraction.h>
#include <src/electronInteraction/nointeraction.h>
#include <src/Hamiltonian/atomichamiltonian.h>
#include <src/Hamiltonian/diatomichamiltonian.h>



VMCApp::VMCApp(const int &myRank, const int &nProcess):
    nProcess(nProcess),
    myRank(myRank)
{
}

//****************************************************************************
void VMCApp::runVMCApp(int nCycles, long idum)
{

    this->nCycles=nCycles/nProcess;
    this->idum = idum - myRank - time(NULL);
    srand(-this->idum);


    setWavefunction();
    setKinetic();
    setPotential();
    setInteraction();
    setHamiltonian();
    setObservables();
    setSolverMethod();
    initilizeAndRunSolver();

}

//****************************************************************************
void VMCApp::setWavefunction(){


    switch (systemType) {
    case ATOMS:
        orbitals = new Hydrogenic(alpha);
        break;

    case MOLECULES:
        orbitals = new Molecular(R, alpha);
        break;
    }

    switch (wavefunctionType) {
    case BASIC:
        jastrow = new NoJastrow;
        break;

    case JASTROW:
        jastrow = new PadeJastrow(nParticles,beta);
        break;
    }

    trialWavefunction = new Wavefunction(cfg, orbitals,jastrow);
}

//****************************************************************************
void VMCApp::setKinetic(){
    kinetic= new Kinetic(trialWavefunction);
}

//****************************************************************************
void VMCApp::setPotential(){
    potential = new CoulombPotential(charge);
}

//****************************************************************************
void VMCApp::setInteraction(){

    switch (InteractionType) {
    case NOINTERACTION:
        electonInteraction = new NoInteraction;
        break;

    case COULOMBINTERACTION:
        electonInteraction = new CoulombInteraction;
        break;
    }
}

//****************************************************************************
void VMCApp::setHamiltonian(){

    switch (systemType) {
    case ATOMS:
        hamiltonian = new AtomicHamiltonian(kinetic,potential,electonInteraction);
        break;

    case MOLECULES:
        hamiltonian = new DiatomicHamiltonian(kinetic,potential,electonInteraction,R);
        break;
    }
}


//****************************************************************************
void VMCApp::setObservables(){
    observables = new Observables(hamiltonian,trialWavefunction);
    observables->loadConfiguration(cfg);
}



//****************************************************************************
void VMCApp::setSolverMethod(){

    switch (solverType) {
    case BF:
        solver = new MCBF(hamiltonian,trialWavefunction,observables);
        break;
    case IS:
        solver =new MCIS(hamiltonian,trialWavefunction,observables);
        break;
    }
}

//****************************************************************************
void VMCApp::initilizeAndRunSolver(){

    solver->loadConfiguration(cfg);
    solver->initializeSolver();
    solver->solve(nCycles,idum);
}

//****************************************************************************
void VMCApp::loadConfiguration(Config *cfg){
    this->cfg=cfg;
    nParticles = cfg->lookup("SolverSettings.N");
    nDimensions = cfg->lookup("SolverSettings.dim");
    charge = cfg->lookup("PotentialSettings.charge");
    minimizationIsEnable =cfg->lookup("setup.minimization");
    blockingIsEnable = cfg->lookup("setup.blocking");
    InteractionType = cfg->lookup("AppSettings.interactionType");
    wavefunctionType = cfg->lookup("AppSettings.wavefunction");
    systemType = cfg->lookup("setup.system");
    solverType= cfg->lookup("AppSettings.solverType");
    R =cfg->lookup("setup.singleRunSettings.R");
    totVariationalDerivate=zeros<vec>(2);
    totEnergyVarDerivate=zeros<vec>(2);
}
//****************************************************************************
double VMCApp::getEnergy(){

    return totEnergy;
}

//****************************************************************************
double VMCApp::getEnergySquared(){

    return totEnergySquared;
}


//****************************************************************************
double VMCApp::getVariance(){

    Variance = totEnergySquared - totEnergy * totEnergy;
    return Variance;
}

//****************************************************************************
double VMCApp::getSigma(){

    Sigma = sqrt(Variance);
    return Sigma;
}

//****************************************************************************
double VMCApp::getAcceptanceRate(){

    return Acceptance;
}

//****************************************************************************
vec VMCApp::getVariationalDerivate(){
    return 2*totEnergyVarDerivate - 2*totVariationalDerivate*totEnergy;
}


//****************************************************************************
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

    if(minimizationIsEnable){
        tmpVec = observables->getVariationalDerivateRatio();
        MPI_Allreduce(&tmpVec[0], &totVariationalDerivate[0], tmpVec.n_rows, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        totVariationalDerivate /= nProcess;

        tmpVec = observables->getEnergyVariationalDerivate();
        MPI_Allreduce(&tmpVec[0], &totEnergyVarDerivate[0], tmpVec.n_rows, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        totEnergyVarDerivate /= nProcess;
    }

    if(blockingIsEnable ){
        observables->writeEnergyVectorToFile(myRank);
    }

    if(myRank==0){
        cout << alpha << ", " << beta << " Energy = " << totEnergy
             << ", Variance = " << getVariance()
             << ", Accepted = " << Acceptance
             << "\n";
    }

}
