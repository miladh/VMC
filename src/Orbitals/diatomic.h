#ifndef DIATOMIC_H
#define DIATOMIC_H
#include<src/Orbitals/orbitals.h>

class Molecular : public Orbitals
{
public:
    Molecular(const double &R, const double &k);
    double orbitalEvaluate(const mat &r, int qNum, int Particle);
    double laplaceOrbitalEvaluate(const mat &r, int qNum, int Particle);
    rowvec gradientOrbitalEvaluate(const mat &r, int qNum, int Particle);
    double getVariationalDerivative(const mat &r, int qNum, int Particle);

    Orbitals* atomicOrbitals;

private:
    double k,R;
    double phi,ddphi;
    mat Rmatrix;
    rowvec dphi;
    double dVariational;
};

#endif // DIATOMIC_H
