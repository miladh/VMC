#include <src/Hamiltonian/atomichamiltonian.h>


AtomicHamiltonian::AtomicHamiltonian()
{
}

/************************************************************
Name:                   getEnergy
Description:            Computes the energy
*/
double AtomicHamiltonian::getEnergy(const mat &r) {

    potentialEnergy = potential->evaluate(r);
    kineticEnergy = kinetic->evaluate(r);
    interactionEnergy = electronInteraction->evaluate(r);
    Energy = interactionEnergy+kineticEnergy+potentialEnergy;

    return Energy;

}
