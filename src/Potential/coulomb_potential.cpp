#include "coulomb_potential.h"

CoulombPotential::CoulombPotential():
    charge(2)
{
}

double CoulombPotential::evaluate(int nParticles,const mat &r){

    double en_energy=electron_nucleus_pot(nParticles,r);
    double ee_energy=electron_electron_pot(nParticles,r);

    return en_energy+ee_energy;

}
double CoulombPotential::electron_nucleus_pot(int nParticles, const mat &r){

    en_potentialEnergy = 0;
    rSingleParticle = 0;

    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = norm(r.row(i),2);
        en_potentialEnergy -= charge / rSingleParticle;
    }

    return en_potentialEnergy;
}


double CoulombPotential::electron_electron_pot(int nParticles, const mat &r){

    rij = 0;
    ee_potentialEnergy = 0;

    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {

            rij= norm( r.row(i)-r.row(j) ,2);
            ee_potentialEnergy += 1 / rij;
        }
    }
    return ee_potentialEnergy;
}
