#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <src/Potential/potential.h>
#include <src/Kinetic/kinetic.h>
#include <src/electronInteraction/electroninteraction.h>

class Hamiltonian
{
public:
    Hamiltonian(Config* cfg, Kinetic* kinetic,Potential* potential,
                ElectronInteraction* electronInteraction);

    virtual vec4 getEnergy(const mat &r)= 0;
    virtual void  setNucleusDistance()=0;

protected:
    Config* cfg;
    Kinetic* kinetic;
    Potential* potential;
    ElectronInteraction* electronInteraction;

    vec4 energyVector;

};

#endif // HAMILTONIAN_H
