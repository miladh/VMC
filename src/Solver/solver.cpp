#include "solver.h"

Solver::Solver(const uint &nParticles, const uint &nDimensions,Hamiltonian *hamiltonian, Wavefunction* TrialWavefunction):
    variationalDerivate(zeros<vec>(2)),
    variationalDerivateSum(zeros<vec>(2)),
    energyVarDerivate(zeros<vec>(2)),
    nParticles(nParticles),
    nDimensions(nDimensions),
    rOld(zeros<mat>(nParticles, nDimensions)),
    rNew(zeros<mat>(nParticles, nDimensions)),
    hamiltonian(hamiltonian),
    TrialWavefunction(TrialWavefunction)
{
}


/************************************************************
Name:               loadConfiguration
Description:        loads different variables
*/
void Solver::loadConfiguration(Config *cfg){
    thermalization=cfg->lookup("AppSettings.thermalization");
    nPreCycles=cfg->lookup("OptimalStepSettings.preCycles");
    minStepLength = cfg->lookup("OptimalStepSettings.minstep");
    maxStepLength = cfg->lookup("OptimalStepSettings.maxstep");
    tolerance = cfg->lookup("OptimalStepSettings.tolerance");
}
