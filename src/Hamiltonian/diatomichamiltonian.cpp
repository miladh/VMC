#include <src/Hamiltonian/diatomichamiltonian.h>

DiatomicHamiltonian::DiatomicHamiltonian(Config* cfg, Kinetic* kinetic,Potential* potential,
                                         ElectronInteraction* electronInteraction, double* R):
    Hamiltonian(cfg, kinetic,potential,electronInteraction),
    R(R)
{
    loadAndSetConfiguration();
}

//*******************************************************************************
vec4 DiatomicHamiltonian::getEnergy(const mat &r) {

    energyVector(1) = potential->evaluate(r-Rmatrix)+potential->evaluate(r+Rmatrix);
    energyVector(2) = kinetic->evaluate(r);
    energyVector(3) = electronInteraction->evaluate(r);
    energyVector(0) = energyVector(1)+energyVector(2)+energyVector(3) + nucleusEnergy;

    return energyVector;

}


//*****************************************************************************
void DiatomicHamiltonian::loadAndSetConfiguration()
{
    nParticles  = cfg->lookup("setup.nParticles");
    nDimensions = cfg->lookup("setup.nDimensions");
    charge      = cfg->lookup("PotentialSettings.charge");

    Rmatrix     = zeros<mat>(nParticles,nDimensions);
    setNucleusDistance();
}

//*****************************************************************************
void DiatomicHamiltonian::setNucleusDistance()
{
    for(int i=0; i < nParticles; i++){
        Rmatrix(i,0) = *R/2;
    }

    nucleusEnergy = charge*charge/(*R);
}
