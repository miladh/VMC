#ifndef BFMINIMIZAER_H
#define BFMINIMIZAER_H
#include <src/Minimizer/minimizer.h>

class BFMinimizer : public Minimizer
{

public:
    BFMinimizer(Config* cfg,const int &myRank, const int &nProcess);
    void runMinimizer();

private:
    double alpha,beta;
    double minAlpha,maxAlpha,minBeta,maxBeta;
    int nVarAlpha,nVarBeta;
    double stepAlpha,stepBeta;

    ofstream myfile;

    void getResultsAndWrite();
    void loadAndSetConfiguration();

};

#endif // BFMINIMIZAER_H
