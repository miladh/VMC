#ifndef BFMINIMIZAER_H
#define BFMINIMIZAER_H
#include <src/Minimizer/minimizer.h>

class BFMinimizer : public Minimizer
{

public:
    BFMinimizer(const int &myRank, const int &nProcess);

    void runMinimizer();
    void loadConfiguration(Config *cfg);

private:
    int nProcess, myRank;

    double alpha,beta;
    double minAlpha,maxAlpha,minBeta,maxBeta;
    int nVarAlpha,nVarBeta;
    double stepAlpha,stepBeta;

    double nCycles;
    long idum;
    double Energy,EnergySquared,Variance, Acceptance,Sigma;

    ofstream myfile;

    void writeToFile();
    VMCApp* vmcapp;
};

#endif // BFMINIMIZAER_H
