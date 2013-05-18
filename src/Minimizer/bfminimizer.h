#ifndef BFMINIMIZAER_H
#define BFMINIMIZAER_H
#include <src/Minimizer/minimizer.h>

class BFMinimizer : public Minimizer
{

public:
    BFMinimizer(Config* cfg,const int &myRank, const int &nProcess);
    void runMinimizer();

private:

    vector<double>alpha,beta,R;
    vec alphaValues,betaValues,RValues;

    ofstream myfile;

    void getResultsAndWrite();
    void loadAndSetConfiguration();



};

#endif // BFMINIMIZAER_H
