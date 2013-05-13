#include <src/Hamiltonian/atomichamiltonian.h>


AtomicHamiltonian::AtomicHamiltonian(Config* cfg, Kinetic* kinetic,Potential* potential,
                                     ElectronInteraction* electronInteraction):
    Hamiltonian(cfg, kinetic,potential,electronInteraction)
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
