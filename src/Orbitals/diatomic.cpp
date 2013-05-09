#include "diatomic.h"
#include <src/Orbitals/hydrogenic.h>

Molecular::Molecular(const double& R, const double& k):
    atomicOrbitals(new Hydrogenic),
    k(k),
    R(R),
    Rmatrix(zeros<mat>(2,3)),
    dphi(zeros<rowvec>(1,3))
{
    Rmatrix(0,0) = R/2.0;
    Rmatrix(1,0) = R/2.0;
    atomicOrbitals->k = k;
}


double Molecular::orbitalEvaluate(const mat &r, int qNum, int Particle){

    phi  = atomicOrbitals->orbitalEvaluate(r-Rmatrix ,qNum, Particle)
         + atomicOrbitals->orbitalEvaluate(r+Rmatrix ,qNum, Particle);
    return phi;

}

double Molecular::laplaceOrbitalEvaluate(const mat &r, int qNum, int Particle){

    ddphi  = atomicOrbitals->laplaceOrbitalEvaluate(r-Rmatrix  ,qNum, Particle)
           + atomicOrbitals->laplaceOrbitalEvaluate(r+Rmatrix  ,qNum, Particle);

    return ddphi;

}

rowvec Molecular::gradientOrbitalEvaluate(const mat &r, int qNum, int Particle){

    dphi  = atomicOrbitals->gradientOrbitalEvaluate(r-Rmatrix ,qNum, Particle)
          + atomicOrbitals->gradientOrbitalEvaluate(r+Rmatrix ,qNum, Particle);

    return dphi;
}

double Molecular::getVariationalDerivative(const mat &, int , int ){
    return 0;
}
