#include <src/Hamiltonian/diatomichamiltonian.h>

DiatomicHamiltonian::DiatomicHamiltonian(Config* cfg, Kinetic* kinetic,Potential* potential,
                                         ElectronInteraction* electronInteraction):
    Hamiltonian(cfg, kinetic,potential,electronInteraction)
{
    loadAndSetConfiguration();
}

//*******************************************************************************
double DiatomicHamiltonian::getEnergy(const mat &r) {

    potentialEnergy = potential->evaluate(r-Rmatrix)+potential->evaluate(r+Rmatrix);
    kineticEnergy   = kinetic->evaluate(r);
    interactionEnergy = electronInteraction->evaluate(r);
    Energy = interactionEnergy + kineticEnergy + potentialEnergy + nucleusEnergy;

    return Energy;

}


//*****************************************************************************
void DiatomicHamiltonian::loadAndSetConfiguration()
{
    nParticles  = cfg->lookup("setup.nParticles");
    nDimensions = cfg->lookup("setup.nDimensions");
    charge      = cfg->lookup("PotentialSettings.charge");
    R           = cfg->lookup("setup.singleRunSettings.R");

    Rmatrix     = zeros<mat>(nParticles,nDimensions);
    for(int i=0; i < nParticles; i++){
        Rmatrix(i,0) = R/2;
    }

    nucleusEnergy = charge*charge/R;

}
