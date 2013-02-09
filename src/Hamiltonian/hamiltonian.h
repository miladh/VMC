#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include"src/Potential/potential.h"
#include"src/Kinetic/kinetic.h"

class Hamiltonian
{
public:
    Hamiltonian();

    double getEnergy(int nDimension, int nParticles, const mat &r);
    Potential* potential;
    Kinetic* kinetic;


protected:
    int nDimension,nParticles,charge;
    double potentialEnergy, kineticEnergy, Energy;


};

#endif // HAMILTONIAN_H
