#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include"src/Potential/potential.h"
#include"src/Kinetic/kinetic.h"

class Hamiltonian
{
public:
    Hamiltonian();

    Potential* potential;
    Kinetic* kinetic;
    double getEnergy(int nParticles, const mat &r);


protected:
    int nDimension,nParticles,charge;
    double potentialEnergy, kineticEnergy, Energy;


};

#endif // HAMILTONIAN_H
