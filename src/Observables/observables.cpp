#include "observables.h"

Observables::Observables(Hamiltonian *hamiltonian, Wavefunction *wavefunction):
    hamiltonian(hamiltonian),
    wavefunction(wavefunction)
{
}


/************************************************************
Name:
Description:
*/
void Observables::calculateObservables()
{

    if(minimize){
        calculateEnergy();
        calculateVariationalDerivateRatio();
    }
    else{
        calculateEnergy();
        if(doBlocking){
            addEnergyToEnergyVector();
        }
    }

}

/************************************************************
Name:
Description:
*/
void Observables::currentConfiguration(const mat& positions)
{
    r = positions;
}

/************************************************************
Name:
Description:
*/
void Observables::calculateEnergy()
{
    deltaE  = hamiltonian->getEnergy(r);
    energy += deltaE;
    energySquared += deltaE*deltaE;
}

/************************************************************
Name:
Description:
*/
double Observables::getEnergy()
{
    return energy/(nCycles-1);
}

/************************************************************
Name:
Description:
*/
double Observables::getEnergySquared()
{
    return energySquared/(nCycles-1);
}

/************************************************************
Name:
Description:
*/
void Observables::calculateVariationalDerivateRatio()
{
    deltaVariationalDerivateRatio = wavefunction->getVariationalDerivate(r);
    for(uint p = 0; p <deltaVariationalDerivateRatio.n_rows; p++ ){
        variationalDerivateRatio(p)  += deltaVariationalDerivateRatio(p);
        energyVariationalDerivate(p) += deltaE*deltaVariationalDerivateRatio(p);
    }
}

/************************************************************
Name:
Description:
*/
vec Observables::getVariationalDerivateRatio()
{
   return variationalDerivateRatio/(nCycles-1);
}

/************************************************************
Name:
Description:
*/
vec Observables::getEnergyVariationalDerivate()
{
   return energyVariationalDerivate/(nCycles-1);

}


/************************************************************
Name:
Description:
*/
void Observables::addEnergyToEnergyVector()
{
    energyVector(cycle) = deltaE;
    cycle += 1;
}

/************************************************************
Name:
Description:
*/
void Observables::writeEnergyVectorToFile(const int& myRank)
{
    ostringstream filename;
    filename << dataPath << dataName << myRank << ".mat";
    energyVector.save(filename.str());
}




/************************************************************
Name:
Description:
*/
void Observables::loadConfiguration(Config *cfg){
    minimize=cfg->lookup("MinimizerSettings.minimize");
    doBlocking= cfg->lookup("BlockingSettings.doBlocking");
    cfg->lookupValue("BlockingSettings.dataPath", dataPath);
    cfg->lookupValue("BlockingSettings.dataName", dataName);
}

/************************************************************
Name:
Description:
*/
void Observables::initializeObservables(const int& nCycles){
    this->nCycles = nCycles;
    cycle=0;
    energy = 0;
    energySquared = 0;
    deltaVariationalDerivateRatio = zeros<vec>(2);
    variationalDerivateRatio = zeros<vec>(2);
    energyVariationalDerivate= zeros<vec>(2);
    energyVector=zeros(nCycles-1);
}
