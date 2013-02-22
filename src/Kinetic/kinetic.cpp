#include "kinetic.h"

Kinetic::Kinetic()
{
}


/************************************************************
Name:               evaluate
Description:        Computes the kinetic energy
*/
double Kinetic::evaluate(const mat &r){

    ddwavefunction = wf->laplace(r);
    KineticEnergy = -0.5 * ddwavefunction;

    return KineticEnergy;
}
