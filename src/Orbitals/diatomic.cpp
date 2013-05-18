#include "diatomic.h"
#include <src/Orbitals/hydrogenic.h>

Diatomic::Diatomic(Config* cfg, Orbitals* orbital, double* R):
    cfg(cfg),
    atomicOrbitals(orbital),
    R(R)
{
    loadAndSetConfiguration();
}

//********************************************************************************
double Diatomic::orbitalEvaluate(const mat &r, int qNum, int Particle)
{
    phi  = atomicOrbitals->orbitalEvaluate(r-Rmatrix ,qNum, Particle)
            + atomicOrbitals->orbitalEvaluate(r+Rmatrix ,qNum, Particle);
    return phi;

}
//********************************************************************************
double Diatomic::laplaceOrbitalEvaluate(const mat &r, int qNum, int Particle)
{
    ddphi  = atomicOrbitals->laplaceOrbitalEvaluate(r-Rmatrix  ,qNum, Particle)
            + atomicOrbitals->laplaceOrbitalEvaluate(r+Rmatrix  ,qNum, Particle);

    return ddphi;

}
//********************************************************************************
rowvec Diatomic::gradientOrbitalEvaluate(const mat &r, int qNum, int Particle)
{
    dphi  = atomicOrbitals->gradientOrbitalEvaluate(r-Rmatrix ,qNum, Particle)
            + atomicOrbitals->gradientOrbitalEvaluate(r+Rmatrix ,qNum, Particle);

    return dphi;
}


//********************************************************************************
double Diatomic::getVariationalDerivative(const mat &r, int qNum, int Particle)
{
    dVariational  = atomicOrbitals->getVariationalDerivative(r-Rmatrix  ,qNum, Particle)
            + atomicOrbitals->getVariationalDerivative(r+Rmatrix  ,qNum, Particle);

    return dVariational;
}


//*****************************************************************************
void  Diatomic::loadAndSetConfiguration()
{
    nParticles  = cfg->lookup("setup.nParticles");
    nDimensions = cfg->lookup("setup.nDimensions");

    dphi = zeros<rowvec>(1,nDimensions);
    Rmatrix     = zeros<mat>(nParticles,nDimensions);
    for(uint i=0; i < nParticles; i++){
        Rmatrix(i,0) = *R/2;
    }


}
