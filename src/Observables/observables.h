#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <src/Hamiltonian/hamiltonian.h>
#include <src/Wavefunction/wavefunction.h>

class Observables
{
public:
    Observables(Config* cfg, Hamiltonian* hamiltonian, Wavefunction* wavefunction);


    void initializeObservables(const int &nCycles);
    void currentConfiguration(const mat& positions);
    void calculateObservables();

    void calculateEnergy();
    void calculateVariationalDerivateRatio();
    void calculateAverageDistance();

    double getEnergy();
    double getEnergySquared();
    double getAverageDistance();
    vec getVariationalDerivateRatio();
    vec getEnergyVariationalDerivate();


    void addEnergyToEnergyVector();
    void writeEnergyVectorToFile(const int &myRank);

private:
    Config* cfg;
    Hamiltonian* hamiltonian;
    Wavefunction* wavefunction;
    int nCycles, cycle;
    int minimize, blockingIsEnable ;
    double deltaE, energy, energySquared;
    double averageDistance;
    vec deltaVariationalDerivateRatio,variationalDerivateRatio,energyVariationalDerivate;
    vec energyVector;
    mat r;
    string dataPath, dataName;


    void loadConfiguration();
};

#endif // OBSERVABLES_H
