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


VMCApp::VMCApp(Config* cfg,const int &myRank, const int &nProcess):
    cfg(cfg),
    nProcess(nProcess),
    myRank(myRank)
{
    loadAndSetConfiguration();
}

//****************************************************************************
void VMCApp::runVMCApp()
{
    setWavefunction();
    setHamiltonian();
    setObservables();
    setAndRunSolver();
    messagePassing();
}

//****************************************************************************
void VMCApp::setWavefunction(){

    switch (systemType) {
    case ATOMS:
        orbitals = new Hydrogenic(alpha);
        break;

    case MOLECULES:
        orbitals = new Molecular(cfg, alpha);
        break;
    }

    switch (wavefunctionType) {
    case BASIC:
        jastrow = new NoJastrow;
        break;

    case JASTROW:
        jastrow = new PadeJastrow(nParticles, beta);
        break;
    }

    trialWavefunction = new Wavefunction(cfg, orbitals,jastrow);
}


//****************************************************************************
void VMCApp::setHamiltonian(){

    Kinetic   *kinetic   = new Kinetic(trialWavefunction);
    Potential *potential = new CoulombPotential(charge);
    ElectronInteraction* electonInteraction;

    switch (InteractionType) {
    case NOINTERACTION:
        electonInteraction = new NoInteraction;
        break;

    case COULOMBINTERACTION:
        electonInteraction = new CoulombInteraction;
        break;
    }


    switch (systemType) {
    case ATOMS:
        hamiltonian = new AtomicHamiltonian(cfg,kinetic,potential,electonInteraction);
        break;

    case MOLECULES:
        hamiltonian = new DiatomicHamiltonian(cfg,kinetic,potential,electonInteraction);
        break;
    }
}


//****************************************************************************
void VMCApp::setObservables(){

    observables = new Observables(cfg,hamiltonian,trialWavefunction);
}


//****************************************************************************
void VMCApp::setAndRunSolver(){

    switch (solverType) {
    case BF:
        solver = new MCBF(cfg, hamiltonian,trialWavefunction,observables);
        break;
    case IS:
        solver = new MCIS(cfg, hamiltonian,trialWavefunction,observables);
        break;
    }

    solver->solve(nCycles,idum);
}


//****************************************************************************
void VMCApp::loadAndSetConfiguration(){

    nParticles  = cfg->lookup("setup.nParticles");
    nDimensions = cfg->lookup("setup.nDimensions");
    charge      = cfg->lookup("PotentialSettings.charge");
    R           = cfg->lookup("setup.singleRunSettings.R");

    minimizationIsEnable = cfg->lookup("setup.minimization");
    blockingIsEnable     = cfg->lookup("setup.blocking");
    InteractionType      = cfg->lookup("AppSettings.interactionType");
    wavefunctionType     = cfg->lookup("AppSettings.wavefunction");

    systemType = cfg->lookup("setup.system");
    solverType = cfg->lookup("AppSettings.solverType");

    nCycles = cfg->lookup("AppSettings.cycles");
    idum    = cfg->lookup("AppSettings.idum");

    nCycles /= nProcess;
    idum -= myRank - time(NULL);
    srand(idum);

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

    tmp = observables->getAverageDistance();
    MPI_Allreduce(&tmp, &averageDistance, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    averageDistance /= nProcess;

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
             << ", Average dist = " << averageDistance
             << "\n";
    }

}
