#include "hamiltonian.h"

Hamiltonian::Hamiltonian()
{
}

/************************************************************
Name:                   getEnergy
Description:            Computes the energy
*/

double Hamiltonian::getEnergy(int nParticles, const mat &r) {

    potentialEnergy = potential->evaluate(nParticles,r);
    kineticEnergy = kinetic->evaluate(nParticles,r);
    Energy = kineticEnergy+potentialEnergy;

    return Energy;

}


