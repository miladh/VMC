#ifndef BFMINIMIZAER_H
#define BFMINIMIZAER_H
#include <src/Minimizer/minimizer.h>

class BFMinimizer : public Minimizer
{

public:
    BFMinimizer(Config* cfg,const int &myRank, const int &nProcess);
    void runMinimizer();

private:
    vec alphaValues,betaValues,RValues;
    int diatomicSystem;
    ofstream myfile;

    void getResultsAndWrite();
    void loadAndSetConfiguration();
    void minimize();



};

#endif // BFMINIMIZAER_H
