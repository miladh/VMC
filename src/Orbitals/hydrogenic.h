#ifndef HYDROGENIC_H
#define HYDROGENIC_H
#include"orbitals.h"

class Hydrogenic : public Orbitals
{
public:
    Hydrogenic();
    double orbitalEvaluate(const mat &r, int qNum, int Particle);
    double LaplaceOrbitalEvaluate(const mat &r, int qNum, int Particle);
    rowvec GradientOrbitalEvaluate(const mat &r, int qNum, int Particle);

private:
    double rNorm;
    double phi,ddphi;
    rowvec dphi;
};

#endif // HYDROGENIC_H
