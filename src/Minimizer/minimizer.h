#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include <src/ConfigurationParser/configurationparser.h>

using namespace arma;
using namespace std;
using namespace libconfig;


class Minimizer
{
public:
    Minimizer(Config* cfg,const int &myRank, const int &nProcess);
    virtual void runMinimizer()=0;

protected:
    Config* cfg;
    uint myRank, nProcess;
    ConfigurationParser* parser;

    double energy,energySquared, variance, acceptance,sigma;
};

#endif // MINIMIZER_H
