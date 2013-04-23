#ifndef CONFIGURATIONPARSER_H
#define CONFIGURATIONPARSER_H

#include <libconfig.h++>
#include <src/Minimizer/minimizer.h>
#include <src/Blocking/blocking.h>
#include <src/OnebodyDensity/onebodydensity.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

using namespace arma;
using namespace std;
using namespace libconfig;



class ConfigurationParser
{
public:
    ConfigurationParser(const int &nProcess, const int &myRank);

    void setup();
    void loadConfiguration();

private:
    int nProcess, myRank;
    int singleRunIsEnable, minimizationIsEnable, blockingIsEnable;
    int minimizerType;
    double alpha, beta;
    double nCycles;
    long idum;
    Config cfg;

    VMCApp* vmcapp;
    Minimizer* minimizer;

    void singleRun();
    void chooseAndRunMinimization();
    void runBlocking();

};

#endif // CONFIGURATIONPARSER_H
