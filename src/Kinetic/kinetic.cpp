#include "kinetic.h"

Kinetic::Kinetic(Config *cfg)
{
    this->cfg=cfg;
}


/************************************************************
Name:               evaluate
Description:        Computes the kinetic energy
*/
double Kinetic::evaluate(int nParticles,const mat &r){

    ddwaveFunction = wf->laplace(nParticles,r,cfg);
    KineticEnergy = -0.5 * ddwaveFunction;

    return KineticEnergy;
}
