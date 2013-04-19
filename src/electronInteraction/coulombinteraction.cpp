#include <src/electronInteraction/coulombinteraction.h>

CoulombInteraction::CoulombInteraction()
{
}


/*********************************************************************************
Name:               electron_electron_pot
Description:        Computes potential energy due to electron-electron interactions
*/

double CoulombInteraction::evaluate(const mat& r){
    rij = 0;
    interactionEnergy = 0;

    for(uint i = 0; i < r.n_rows; i++) {
        for(uint j = i + 1; j < r.n_rows; j++) {
            rij= norm( r.row(i)-r.row(j) ,2);
            interactionEnergy += 1 / rij;
        }
    }
    return interactionEnergy;
}
