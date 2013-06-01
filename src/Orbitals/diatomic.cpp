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
    phi  = atomicOrbitals->orbitalEvaluate(r+Rmatrix ,qNum/2, Particle)+
            signFunc(qNum)*atomicOrbitals->orbitalEvaluate(r-Rmatrix ,qNum/2, Particle);
    return phi;
}
//********************************************************************************
double Diatomic::laplaceOrbitalEvaluate(const mat &r, int qNum, int Particle)
{
    ddphi  = atomicOrbitals->laplaceOrbitalEvaluate(r+Rmatrix  ,qNum/2, Particle)+
            signFunc(qNum)* atomicOrbitals->laplaceOrbitalEvaluate(r-Rmatrix  ,qNum/2, Particle);

    return ddphi;

}
//********************************************************************************
rowvec Diatomic::gradientOrbitalEvaluate(const mat &r, int qNum, int Particle)
{
    dphi  = atomicOrbitals->gradientOrbitalEvaluate(r+Rmatrix ,qNum/2, Particle)+
            signFunc(qNum)* atomicOrbitals->gradientOrbitalEvaluate(r-Rmatrix ,qNum/2, Particle);

    return dphi;
}


//********************************************************************************
double Diatomic::getVariationalDerivative(const mat &r, int qNum, int Particle)
{
    dVariational  = atomicOrbitals->getVariationalDerivative(r+Rmatrix  ,qNum/2, Particle)+
            signFunc(qNum)* atomicOrbitals->getVariationalDerivative(r-Rmatrix  ,qNum/2, Particle);

    return dVariational;
}


//*****************************************************************************
void  Diatomic::setNucleusDistance()
{
    for(uint i=0; i < nParticles; i++){
        Rmatrix(i,0) = *R/2;
    }
}


//*****************************************************************************
void  Diatomic::loadAndSetConfiguration()
{
    nParticles  = cfg->lookup("setup.nParticles");
    nDimensions = cfg->lookup("setup.nDimensions");

    dphi = zeros<rowvec>(1,nDimensions);
    Rmatrix     = zeros<mat>(nParticles,nDimensions);

    setNucleusDistance();
}


