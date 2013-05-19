#include "observables.h"
#include <src/Hamiltonian/atomichamiltonian.h>
#include <src/Hamiltonian/diatomichamiltonian.h>

Observables::Observables(Config* cfg,Hamiltonian *hamiltonian, Wavefunction *wavefunction):
    cfg(cfg),
    hamiltonian(hamiltonian),
    wavefunction(wavefunction)
{
    loadConfiguration();
}


//************************************************************
void Observables::calculateObservables()
{

    if(minimize){
        calculateEnergy();
        calculateVariationalDerivateRatio();
    }
    else{
        calculateEnergy();
        calculateAverageDistance();
        addPositionsToPositionMatrix();

        if(blockingIsEnable){
            addEnergyToEnergyVector();
        }
    }

}

//************************************************************
void Observables::initializeObservables(const int& nCycles){
    this->nCycles = nCycles;
    cycle  = 0;
    energy = 0;
    energySquared = 0;
    averageDistance = 0;
    deltaVariationalDerivateRatio = zeros<vec>(2);
    variationalDerivateRatio  = zeros<vec>(2);
    energyVariationalDerivate = zeros<vec>(2);
    energyVector = zeros(nCycles-1);
}

//************************************************************
void Observables::currentConfiguration(const mat& positions)
{
    r = positions;
}


//************************************************************
void Observables::calculateEnergy()
{
    deltaE  = hamiltonian->getEnergy(r);
    energy += deltaE;
    energySquared += deltaE*deltaE;
}

//************************************************************
void Observables::calculateAverageDistance()
{
    for (uint i=0; i<r.n_rows; i++) {
        for (uint j=i+1; j<r.n_rows; j++) {
            averageDistance+= norm( r.row(i)-r.row(j) ,2);
        }
    }
}
//************************************************************
void Observables::calculateVariationalDerivateRatio()
{
    deltaVariationalDerivateRatio = wavefunction->getVariationalDerivate(r);
    for(uint p = 0; p <deltaVariationalDerivateRatio.n_rows; p++ ){
        variationalDerivateRatio(p)  += deltaVariationalDerivateRatio(p);
        energyVariationalDerivate(p) += deltaE*deltaVariationalDerivateRatio(p);
    }
}

//************************************************************6
double Observables::getEnergy()
{
    return energy/(nCycles-1);
}

//************************************************************
double Observables::getEnergySquared()
{
    return energySquared/(nCycles-1);
}

//************************************************************
double Observables::getAverageDistance()
{
    return averageDistance/(nCycles-1);
}
//************************************************************
vec Observables::getVariationalDerivateRatio()
{
    return variationalDerivateRatio/(nCycles-1);
}


//************************************************************
vec Observables::getEnergyVariationalDerivate()
{
    return energyVariationalDerivate/(nCycles-1);

}

//************************************************************
void Observables::addPositionsToPositionMatrix()
{
    positionsMat.push_back(r);
}

//************************************************************
void Observables::addEnergyToEnergyVector()
{
    energyVector(cycle) = deltaE;
    cycle += 1;
}


//************************************************************
void Observables::writePositionMatrixToFile(const int& myRank)
{
    ofstream myfile;
    ostringstream filename;
    filename << "../vmc/DATA/onebodyDensity/OBD"<< myRank << ".mat";
    myfile.open(filename.str().c_str());
    for(uint i=0; i<positionsMat.size(); i++){
        myfile << positionsMat[i] << endl;
    }
}

//************************************************************
void Observables::writeEnergyVectorToFile(const int& myRank)
{
    ostringstream filename;
    filename << dataPath << dataName << myRank << ".mat";
    energyVector.save(filename.str());
}



//************************************************************
void Observables::loadConfiguration(){

    minimize         = cfg->lookup("setup.minimization");
    blockingIsEnable = cfg->lookup("setup.blocking");
    cfg->lookupValue("setup.BlockingSettings.dataPath", dataPath);
    cfg->lookupValue("setup.BlockingSettings.dataName", dataName);
}
