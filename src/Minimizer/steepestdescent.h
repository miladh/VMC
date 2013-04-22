#ifndef STEEPESTDESCENT_H
#define STEEPESTDESCENT_H

#include <src/Minimizer/minimizer.h>


class SteepestDescent : public Minimizer
{
public:
    SteepestDescent(const int &myRank, const int &nProcess);


    void runMinimizer();
    void loadConfiguration(Config *cfg);

private:
    int nProcess, myRank;

    double alpha,beta;
    double stepAlpha,stepBeta;

    double nCycles;
    long idum;
    double energy,energySquared,variance, acceptance,sigma;

    VMCApp* vmcapp;
    vec variationalDerivate;

};

#endif // STEEPESTDESCENT_H
