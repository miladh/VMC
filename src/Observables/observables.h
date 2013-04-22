#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <src/Hamiltonian/hamiltonian.h>
#include <src/Wavefunction/wavefunction.h>

class Observables
{
public:
    Observables(Hamiltonian* hamiltonian, Wavefunction* wavefunction);

    void loadConfiguration(Config *cfg);
    void initializeObservables(const int &nCycles);
    void calculateObservables();

    void currentConfiguration(const mat& positions);

    void calculateEnergy();
    double getEnergy();
    double getEnergySquared();

    void calculateVariationalDerivateRatio();
    vec getVariationalDerivateRatio();
    vec getEnergyVariationalDerivate();

    void addEnergyToEnergyVector();
    void writeEnergyVectorToFile(const int &myRank);

private:
    Hamiltonian* hamiltonian;
    Wavefunction* wavefunction;
    int nCycles, cycle;
    int minimize, doBlocking;
    double deltaE, energy, energySquared;
    vec deltaVariationalDerivateRatio,variationalDerivateRatio,energyVariationalDerivate;
    vec energyVector;
    mat r;
    string dataPath, dataName;

};

#endif // OBSERVABLES_H
