#include "diatomic.h"
#include <src/Orbitals/hydrogenic.h>

Molecular::Molecular(Config* cfg, const double& k):
    cfg(cfg),
    atomicOrbitals(new Hydrogenic(k))
{
    loadAndSetConfiguration();
}

//********************************************************************************
double Molecular::orbitalEvaluate(const mat &r, int qNum, int Particle)
{
    phi  = atomicOrbitals->orbitalEvaluate(r-Rmatrix ,qNum, Particle)
            + atomicOrbitals->orbitalEvaluate(r+Rmatrix ,qNum, Particle);
    return phi;

}
//********************************************************************************
double Molecular::laplaceOrbitalEvaluate(const mat &r, int qNum, int Particle)
{
    ddphi  = atomicOrbitals->laplaceOrbitalEvaluate(r-Rmatrix  ,qNum, Particle)
            + atomicOrbitals->laplaceOrbitalEvaluate(r+Rmatrix  ,qNum, Particle);

    return ddphi;

}
//********************************************************************************
rowvec Molecular::gradientOrbitalEvaluate(const mat &r, int qNum, int Particle)
{
    dphi  = atomicOrbitals->gradientOrbitalEvaluate(r-Rmatrix ,qNum, Particle)
            + atomicOrbitals->gradientOrbitalEvaluate(r+Rmatrix ,qNum, Particle);

    return dphi;
}


//********************************************************************************
double Molecular::getVariationalDerivative(const mat &r, int qNum, int Particle)
{
    dVariational  = atomicOrbitals->getVariationalDerivative(r-Rmatrix  ,qNum, Particle)
            + atomicOrbitals->getVariationalDerivative(r+Rmatrix  ,qNum, Particle);

    return dVariational;
}


//*****************************************************************************
void  Molecular::loadAndSetConfiguration()
{
    nParticles  = cfg->lookup("setup.nParticles");
    nDimensions = cfg->lookup("setup.nDimensions");
    R           = cfg->lookup("setup.singleRunSettings.R");

    dphi = zeros<rowvec>(1,nDimensions);
    Rmatrix     = zeros<mat>(nParticles,nDimensions);
    for(uint i=0; i < nParticles; i++){
        Rmatrix(i,0) = R/2;
    }


}
