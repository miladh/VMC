#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include "src/VMCApp/vmcapp.h"

using namespace arma;
using namespace std;
using namespace libconfig;


class Minimizer
{
public:
    Minimizer(const int &myRank, const int &nProcess);

    void runMinimizaer();
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

#endif // MINIMIZER_H
