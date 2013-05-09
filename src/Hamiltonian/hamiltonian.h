#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <src/Potential/potential.h>
#include <src/Kinetic/kinetic.h>
#include <src/electronInteraction/electroninteraction.h>

class Hamiltonian
{
public:
    Hamiltonian(Kinetic* kinetic,Potential* potential,
                ElectronInteraction* electronInteraction);

    virtual double getEnergy(const mat &r)= 0;

protected:
    Kinetic* kinetic;
    Potential* potential;
    ElectronInteraction* electronInteraction;

};

#endif // HAMILTONIAN_H
