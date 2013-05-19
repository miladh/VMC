#ifndef DIATOMIC_H
#define DIATOMIC_H
#include<src/Orbitals/orbitals.h>

class Diatomic : public Orbitals
{
public:
    Diatomic(Config *cfg, Orbitals *orbital, double *R);
    double orbitalEvaluate(const mat &r, int qNum, int Particle);
    double laplaceOrbitalEvaluate(const mat &r, int qNum, int Particle);
    double getVariationalDerivative(const mat &r, int qNum, int Particle);
    rowvec gradientOrbitalEvaluate(const mat &r, int qNum, int Particle);

    void  setNucleusDistance();

private:
    Config* cfg;
    Orbitals* atomicOrbitals;
    uint nParticles, nDimensions;
    double* R;
    double phi,ddphi,dVariational;
    mat Rmatrix;
    rowvec dphi;


    void loadAndSetConfiguration();
};

#endif // DIATOMIC_H
