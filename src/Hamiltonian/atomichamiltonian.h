#ifndef ATOMICHAMILTONIAN_H
#define ATOMICHAMILTONIAN_H

#include <src/Hamiltonian/hamiltonian.h>


class AtomicHamiltonian : public Hamiltonian
{
public:
    AtomicHamiltonian(Config *cfg, Kinetic *kinetic, Potential *potential,
                      ElectronInteraction *electronInteraction);

    double getEnergy(const mat &r);

private:
    double potentialEnergy, kineticEnergy,interactionEnergy, Energy;

};

#endif // ATOMICHAMILTONIAN_H
