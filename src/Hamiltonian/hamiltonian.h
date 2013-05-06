#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <src/Potential/potential.h>
#include <src/Kinetic/kinetic.h>
#include <src/electronInteraction/electroninteraction.h>

class Hamiltonian
{
public:
    Hamiltonian();

    virtual double getEnergy(const mat &r)= 0;

    Potential* potential;
    Kinetic* kinetic;
    ElectronInteraction* electronInteraction;

};

#endif // HAMILTONIAN_H
