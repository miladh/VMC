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
    void writeToFile();
    Config *cfg;
    VMCApp* vmcapp;
    ofstream myfile;



private:
    double alpha,beta;
    double minAlpha,maxAlpha,minBeta,maxBeta;
    int nVarAlpha,nVarBeta;
    double stepAlpha,stepBeta;

    int nCycles;
    long idum;
    int nProcess, myRank;
    double Energy,EnergySquared,Variance, Acceptance,Sigma;

};

#endif // MINIMIZER_H
