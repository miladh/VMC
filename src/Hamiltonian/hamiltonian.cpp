#include "hamiltonian.h"
#include <src/electronInteraction/coulombinteraction.h>
#include <src/electronInteraction/nointeraction.h>

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
    interactionEnergy = electronInteraction->evaluate(r);
    Energy = interactionEnergy+kineticEnergy+potentialEnergy;

    return Energy;

}


