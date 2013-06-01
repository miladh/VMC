#include "coulombPotential.h"

CoulombPotential::CoulombPotential(const int &charge):
    charge(charge)
{
}

//******************************************************************************
double CoulombPotential::evaluate(const mat &r){
    return electronNucleusPotential(r);
}

//********************************************************************************
double CoulombPotential::electronNucleusPotential(const mat &r){
    enPotentialEnergy = 0;

    for(uint i = 0; i < r.n_rows; i++) {
        rSingleParticle = norm(r.row(i),2);
        enPotentialEnergy -= charge / rSingleParticle;
    }

    return enPotentialEnergy;
}

