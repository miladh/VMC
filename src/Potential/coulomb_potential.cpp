#include "coulomb_potential.h"

CoulombPotential::CoulombPotential():
    charge(2)
{
}

double CoulombPotential::evaluate(int nDimensions, int nParticles,const mat &r){

    double en_energy=electron_nucleus_pot(nDimensions,nParticles,r);
    double ee_energy=electron_electron_pot(nDimensions,nParticles,r);

    return en_energy+ee_energy;

}
double CoulombPotential::electron_nucleus_pot(int nDimensions, int nParticles,const mat &r){

    double en_potentialEnergy = 0;
    double rSingleParticle = 0;

    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        en_potentialEnergy -= charge / sqrt(rSingleParticle);
    }

    return en_potentialEnergy;
}


double CoulombPotential::electron_electron_pot(int nDimensions, int nParticles,const mat &r){

    double r12 = 0;
    double ee_potentialEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            ee_potentialEnergy += 1 / sqrt(r12);
        }
    }
    return ee_potentialEnergy;
}
