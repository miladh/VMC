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
    void writeToFile();
    double McIntegrator();
    void normalize();
    void setWavefunction();

private:
    Wavefunction* wf;
    uint nDimensions,nParticles;
    int charge;
    double alpha,beta;
    double nSteps,dr, a, b;
    double nCycles,wfValue;

    int nNodes, myRank;
    int WavefunctionType;
    long idum;
    ofstream myfile;

    mat r;
    vec rho;
};

#endif // ONEBODYDENSITY_H
