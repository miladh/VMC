#include <src/Hamiltonian/atomichamiltonian.h>


AtomicHamiltonian::AtomicHamiltonian(Config* cfg, Kinetic* kinetic,Potential* potential,
                                     ElectronInteraction* electronInteraction):
    Hamiltonian(cfg, kinetic,potential,electronInteraction)
{
}

//***********************************************************************************
vec4 AtomicHamiltonian::getEnergy(const mat &r) {

    energyVector(1) = potential->evaluate(r);
    energyVector(2) = kinetic->evaluate(r);
    energyVector(3) = electronInteraction->evaluate(r);
    energyVector(0) = energyVector(1)+energyVector(2)+energyVector(3);

    return energyVector;

}
