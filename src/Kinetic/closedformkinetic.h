#ifndef CLOSEDFORMKINETIC_H
#define CLOSEDFORMKINETIC_H

#include "src/Kinetic/kinetic.h"

class ClosedFormKinetic : public Kinetic
{
public:
    ClosedFormKinetic();
    double evaluate(int nParticles, const mat &r);
    double localEnergy;


private:
    int charge;
    double r1, r2, rij,r1r2;
    double E_L1,eIntEnergy,eContributor;
};

#endif // CLOSEDFORMKINETIC_H
