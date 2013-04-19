#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include"src/Potential/potential.h"
#include"src/Kinetic/kinetic.h"
#include <src/electronInteraction/electroninteraction.h>

class Hamiltonian
{
public:
    Hamiltonian();

    double getEnergy(const mat &r);
    Potential* potential;
    Kinetic* kinetic;
    ElectronInteraction* electronInteraction;

private:
    double potentialEnergy, kineticEnergy,interactionEnergy, Energy;

};

#endif // HAMILTONIAN_H
