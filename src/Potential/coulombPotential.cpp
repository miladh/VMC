#include "coulombPotential.h"

CoulombPotential::CoulombPotential()
{
}

/******************************************************************************
Name:               evaluate
Description:        Computes the total potential energy
*/

double CoulombPotential::evaluate(const mat &r){
    return electronNucleusPotential(r);
}

/********************************************************************************
Name:               electron_nucleus_pot
Description:        Computes potential energy due to electron-nucleus interactions
*/


double CoulombPotential::electronNucleusPotential(const mat &r){
    enPotentialEnergy = 0;

    for(uint i = 0; i < r.n_rows; i++) {
        rSingleParticle = norm(r.row(i),2);
        enPotentialEnergy -= charge / rSingleParticle;
    }

    return enPotentialEnergy;
}

