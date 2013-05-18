#ifndef DIATOMICHAMILTONIAN_H
#define DIATOMICHAMILTONIAN_H

#include <src/Hamiltonian/hamiltonian.h>

class DiatomicHamiltonian : public Hamiltonian
{
public:
    DiatomicHamiltonian(Config *cfg, Kinetic* kinetic, Potential* potential,
                        ElectronInteraction* electronInteraction,double* R);

    double getEnergy(const mat &r);

private:
    double* R;
    mat Rmatrix;
    double potentialEnergy, kineticEnergy,interactionEnergy, nucleusEnergy,Energy;
    int nParticles,nDimensions,charge;

    void loadAndSetConfiguration();
};

#endif // DIATOMICHAMILTONIAN_H
