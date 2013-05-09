#include <src/Hamiltonian/diatomichamiltonian.h>

DiatomicHamiltonian::DiatomicHamiltonian(const double& R):
    R(R),
    Rmatrix(zeros<mat>(2,3))
{
    Rmatrix(0,0) = R/2;
    Rmatrix(1,0) = R/2;

}

/************************************************************
Name:                   getEnergy
Description:            Computes the energy
*/

double DiatomicHamiltonian::getEnergy(const mat &r) {

    potentialEnergy = potential->evaluate(r-Rmatrix)+potential->evaluate(r+Rmatrix);
    kineticEnergy   = kinetic->evaluate(r);
    interactionEnergy = electronInteraction->evaluate(r);
    Energy = interactionEnergy+kineticEnergy+potentialEnergy + 1/R;

    return Energy;

}
