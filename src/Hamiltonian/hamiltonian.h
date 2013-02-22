#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include"src/Potential/potential.h"
#include"src/Kinetic/kinetic.h"

class Hamiltonian
{
public:
    Hamiltonian();

    double getEnergy(const mat &r);
    Potential* potential;
    Kinetic* kinetic;

private:
    double potentialEnergy, kineticEnergy, Energy;

};

#endif // HAMILTONIAN_H
