#include "coulomb_potential.h"

CoulombPotential::CoulombPotential()
{
}

/******************************************************************************
Name:               evaluate
Description:        Computes the total potential energy
*/

double CoulombPotential::evaluate(const mat &r){
    return electronNucleusPotential(r)+electronElectronPotential(r);
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


/*********************************************************************************
Name:               electron_electron_pot
Description:        Computes potential energy due to electron-electron interactions
*/

double CoulombPotential::electronElectronPotential(const mat &r){
    rij = 0;
    eePotentialEnergy = 0;

    for(uint i = 0; i < r.n_rows; i++) {
        for(uint j = i + 1; j < r.n_rows; j++) {

            rij= norm( r.row(i)-r.row(j) ,2);
            eePotentialEnergy += 1 / rij;
        }
    }
    return eePotentialEnergy;
}
