#include "solver.h"

Solver::Solver()
{
}


/************************************************************
Name:               loadConfiguration
Description:        loads different variables
*/
void Solver::loadConfiguration(Config *cfg){
    nDimensions=cfg->lookup("SolverSettings.dim");
    nParticles=cfg->lookup("SolverSettings.N");
    thermalization=cfg->lookup("AppSettings.thermalization");

    nPreCycles=cfg->lookup("OptimalStepSettings.preCycles");
    minStepLength = cfg->lookup("OptimalStepSettings.minstep");
    maxStepLength = cfg->lookup("OptimalStepSettings.maxstep");
    tolerance = cfg->lookup("OptimalStepSettings.tolerance");


}
