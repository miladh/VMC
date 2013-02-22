#ifndef COULOMB_POTENTIAL_H
#define COULOMB_POTENTIAL_H

#include "src/Potential/potential.h"

class CoulombPotential : public Potential
{
public:
    CoulombPotential();
    double evaluate(const mat &r);

private:
    double electronNucleusPotential(const mat &r);
    double electronElectronPotential(const mat &r);

    double enPotentialEnergy,eePotentialEnergy;
    double rSingleParticle,rij;
};

#endif // COULOMB_POTENTIAL_H
