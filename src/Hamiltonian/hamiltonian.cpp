#include "hamiltonian.h"

Hamiltonian::Hamiltonian()
{
}

/************************************************************
Name:                   getEnergy
Description:            Computes the energy
*/

double Hamiltonian::getEnergy(const mat &r) {

    potentialEnergy = potential->evaluate(r);
    kineticEnergy = kinetic->evaluate(r);
    Energy = kineticEnergy+potentialEnergy;

    return Energy;

}


