#ifndef HYDROGENIC_H
#define HYDROGENIC_H
#include<src/Orbitals/orbitals.h>

class Hydrogenic : public Orbitals
{
public:
    Hydrogenic(double *k);
    double orbitalEvaluate(const mat &r, int qNum, int Particle);
    double laplaceOrbitalEvaluate(const mat &r, int qNum, int Particle);
    double getVariationalDerivative(const mat &r, int qNum, int Particle);
    rowvec gradientOrbitalEvaluate(const mat &r, int qNum, int Particle);

    void setNucleusDistance(){}

private:
    double* k;
    double rNorm;
    double phi,ddphi;
    rowvec dphi;
    double dVariational;
};

#endif // HYDROGENIC_H
