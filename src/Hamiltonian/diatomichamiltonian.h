#ifndef DIATOMICHAMILTONIAN_H
#define DIATOMICHAMILTONIAN_H

#include <src/Hamiltonian/hamiltonian.h>

class DiatomicHamiltonian : public Hamiltonian
{
public:
    DiatomicHamiltonian(Config *cfg, Kinetic* kinetic, Potential* potential,
                        ElectronInteraction* electronInteraction,double* R);

    vec4 getEnergy(const mat &r);
    void  setNucleusDistance();

private:
    double* R;
    mat Rmatrix;
    double nucleusEnergy;
    int nParticles,nDimensions,charge;

    void loadAndSetConfiguration();
};

#endif // DIATOMICHAMILTONIAN_H
