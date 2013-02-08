#ifndef COULOMB_POTENTIAL_H
#define COULOMB_POTENTIAL_H

#include "src/Potential/potential.h"

class CoulombPotential : public Potential
{
public:
    CoulombPotential();
    double evaluate(int nDimensions, int nParticles,const mat &r);

private:
    double electron_nucleus_pot(int nDimensions, int nParticles,const mat &r);
    double electron_electron_pot(int nDimensions, int nParticles,const mat &r);
    int charge;

};

#endif // COULOMB_POTENTIAL_H
