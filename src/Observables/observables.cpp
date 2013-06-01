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
        if(cycle %5==0){
            addPositionsToPositionMatrix();
        }

        if(blockingIsEnable){
            addEnergyTototEnergyVector();
        }
    }

}

//************************************************************
void Observables::initializeObservables(const int& nCycles){
    this->nCycles = nCycles;
    cycle  = 0;
    energySquared = 0;
    averageDistance = 0;
    deltaVariationalDerivateRatio = zeros<vec>(2);
    variationalDerivateRatio  = zeros<vec>(2);
    energyVariationalDerivate = zeros<vec>(2);
    totEnergyVector = zeros(nCycles-1);
    energyVector = zeros(4);
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
    energyVector += deltaE;
    energySquared += deltaE(0)*deltaE(0);
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
        energyVariationalDerivate(p) += deltaE(0)*deltaVariationalDerivateRatio(p);
    }
}

//************************************************************6
vec4 Observables::getEnergy()
{
    return energyVector/(nCycles-1);
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
void Observables::addEnergyTototEnergyVector()
{
    totEnergyVector(cycle) = deltaE(0);
    cycle += 1;
}


//************************************************************
void Observables::writePositionMatrixToFile(const int& myRank)
{
    ostringstream filename;
    filename << "../vmc/DATA/onebodyDensity/OBD"<< myRank <<".bin";

    ofstream myfile (filename.str().c_str(), ios::out | ios::binary);

    for(uint i=0; i<positionsMat.size(); i++){
        for(uint j=0; j < positionsMat[i].n_rows; j++){
            for(uint k=0; k<positionsMat[i].n_cols; k++){

                myfile.write((char*)&positionsMat[i](j,k), sizeof(double));
            }
        }
    }
}

//************************************************************
void Observables::writetotEnergyVectorToFile(const int& myRank)
{
    ostringstream filename;
    filename << dataPath << dataName << myRank << ".mat";
    totEnergyVector.save(filename.str());
}



//************************************************************
void Observables::loadConfiguration(){

    minimize         = cfg->lookup("setup.minimization");
    blockingIsEnable = cfg->lookup("setup.blocking");
    cfg->lookupValue("setup.BlockingSettings.dataPath", dataPath);
    cfg->lookupValue("setup.BlockingSettings.dataName", dataName);
}
