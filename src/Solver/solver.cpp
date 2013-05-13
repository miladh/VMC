#include "solver.h"

Solver::Solver(Config* cfg,Hamiltonian* hamiltonian, Wavefunction* trialWavefunction, Observables* observables):
    cfg(cfg),
    hamiltonian(hamiltonian),
    trialWavefunction(trialWavefunction),
    observables(observables)
{
    loadAndSetConfiguration();
}


/************************************************************
Name:               loadConfiguration
Description:        loads different variables
*/
void Solver::loadAndSetConfiguration(){

    nParticles      = cfg->lookup("setup.nParticles");
    nDimensions     = cfg->lookup("setup.nDimensions");

    D               = cfg->lookup("solverSettings.IS.D");
    timeStep        = cfg->lookup("solverSettings.IS.timeStep");

    thermalization  = cfg->lookup("AppSettings.thermalization");
    nPreCycles      = cfg->lookup("solverSettings.BF.OptimalStepSettings.preCycles");
    minStepLength   = cfg->lookup("solverSettings.BF.OptimalStepSettings.minstep");
    maxStepLength   = cfg->lookup("solverSettings.BF.OptimalStepSettings.maxstep");
    tolerance       = cfg->lookup("solverSettings.BF.OptimalStepSettings.tolerance");

    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
}
