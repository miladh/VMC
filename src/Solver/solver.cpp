#include "solver.h"

Solver::Solver(Hamiltonian *hamiltonian, Wavefunction* trialWavefunction, Observables* observables):
    hamiltonian(hamiltonian),
    trialWavefunction(trialWavefunction),
    observables(observables)
{
}


/************************************************************
Name:               loadConfiguration
Description:        loads different variables
*/
void Solver::loadConfiguration(Config *cfg){
    nParticles = cfg->lookup("SolverSettings.N");
    nDimensions = cfg->lookup("SolverSettings.dim");
    thermalization=cfg->lookup("AppSettings.thermalization");
    nPreCycles=cfg->lookup("OptimalStepSettings.preCycles");
    minStepLength = cfg->lookup("OptimalStepSettings.minstep");
    maxStepLength = cfg->lookup("OptimalStepSettings.maxstep");
    tolerance = cfg->lookup("OptimalStepSettings.tolerance");
}


/************************************************************
Name:               loadConfiguration
Description:        loads different variables
*/
void Solver::initializeSolver(){
    rOld=zeros<mat>(nParticles, nDimensions);
    rNew=zeros<mat>(nParticles, nDimensions);

}
