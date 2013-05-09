#ifndef DIATOMIC_H
#define DIATOMIC_H
#include<src/Orbitals/orbitals.h>

class Molecular : public Orbitals
{
public:
    Molecular(const double &R, const double &k);
    double orbitalEvaluate(const mat &r, int qNum, int Particle);
    double laplaceOrbitalEvaluate(const mat &r, int qNum, int Particle);
    double getVariationalDerivative(const mat &r, int qNum, int Particle);
    rowvec gradientOrbitalEvaluate(const mat &r, int qNum, int Particle);


private:
    Orbitals* atomicOrbitals;
    double R;
    double phi,ddphi,dVariational;
    mat Rmatrix;
    rowvec dphi;
};

#endif // DIATOMIC_H
