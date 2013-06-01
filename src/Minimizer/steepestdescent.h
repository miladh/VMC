#ifndef STEEPESTDESCENT_H
#define STEEPESTDESCENT_H

#include <src/Minimizer/minimizer.h>


class SteepestDescent : public Minimizer
{
public:
    SteepestDescent(Config* cfg, const int &myRank, const int &nProcess);
    void runMinimizer();


private:
    double step, epsilon;
    int nVariationalParameters, maxIteration;
    int diatomicSystem;

    ofstream myfile;
    vec variationalDerivate;

    void loadAndSetConfiguration();
    void  minimize();
    void getResultsAndWrite();
    int signFunc(double varDer);

};

#endif // STEEPESTDESCENT_H
