#ifndef DIATOMICHAMILTONIAN_H
#define DIATOMICHAMILTONIAN_H

#include <src/Hamiltonian/hamiltonian.h>

class DiatomicHamiltonian : public Hamiltonian
{
public:
    DiatomicHamiltonian(const double &R);

    double getEnergy(const mat &r);

private:
    double R;
    mat Rmatrix;
    double potentialEnergy, kineticEnergy,interactionEnergy, Energy;
    uint charge;
};

#endif // DIATOMICHAMILTONIAN_H
