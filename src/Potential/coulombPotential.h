#ifndef COULOMB_POTENTIAL_H
#define COULOMB_POTENTIAL_H

#include <src/Potential/potential.h>

class CoulombPotential : public Potential
{
public:
    CoulombPotential(const int &charge);
    double evaluate(const mat &r);

private:
    double electronNucleusPotential(const mat &r);
    double enPotentialEnergy,eePotentialEnergy;
    double rSingleParticle,rij;
    int charge;
};

#endif // COULOMB_POTENTIAL_H
