#include "hamiltonian.h"

Hamiltonian::Hamiltonian()
{
}

double Hamiltonian::getEnergy(int nDimension,int nParticles, const mat &r) {

    potentialEnergy = potential->evaluate(nDimension,nParticles,r);
    kineticEnergy = kinetic->evaluate(nDimension,nParticles,r);
    Energy = kineticEnergy;
//    Energy = kineticEnergy+potentialEnergy;

    return Energy;

}


