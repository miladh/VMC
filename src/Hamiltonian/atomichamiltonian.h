#ifndef ATOMICHAMILTONIAN_H
#define ATOMICHAMILTONIAN_H

#include <src/Hamiltonian/hamiltonian.h>


class AtomicHamiltonian : public Hamiltonian
{
public:
    AtomicHamiltonian(Config *cfg, Kinetic *kinetic, Potential *potential,
                      ElectronInteraction *electronInteraction);

    vec4 getEnergy(const mat &r);
    void setNucleusDistance(){}

};

#endif // ATOMICHAMILTONIAN_H
