#ifndef VMCAPP_H
#define VMCAPP_H

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



class VMCApp
{
public:
    VMCApp(const int &nProcess, const int &myRank);
    void options();

private:
    int nProcess, myRank;
    Config cfg;

    void loadConfiguration();
    void singleRun();
    void chooseAndRunMinimization();
    void runBlocking();

};

#endif // VMCAPP_H
