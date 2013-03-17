#ifndef ONEBODYDENSITY_H
#define ONEBODYDENSITY_H


#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include "src/Jastrow/jastrow.h"
#include "src/slater/slater.h"
#include "src/Wavefunction/wavefunction.h"


using namespace arma;
using namespace std;
using namespace libconfig;

class OnebodyDensity
{
public:
    OnebodyDensity(const int &nProcess, const int &myRank);

    void computeOnebodyDensity();
    void loadConfiguration(Config *cfg);
    void writeToFile(mat r, vec rho);
    double McIntegrator(mat r);
    vec normalize(vec rho);
    Wavefunction* setWavefunction();

private:
    double nCycles;
    Wavefunction* wf;
    uint nDimensions;
    uint nParticles;
    long idum;
    int nNodes, myRank;

    double alpha,beta;
    double nSteps,dr, a, b;
    double wfValue;
    ofstream myfile;
    int charge;
    int WavefunctionType;
};

#endif // ONEBODYDENSITY_H
