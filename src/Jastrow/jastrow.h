#ifndef JASTROW_H
#define JASTROW_H


#include <armadillo>
#include <iostream>
#include <libconfig.h++>


using namespace arma;
using namespace std;
using namespace libconfig;

class Jastrow
{
public:
    Jastrow();
    void setaValues(const uint &nParticles);
    double evaluateJastrow(const mat &r);
    double laplaceJastrowEvaluate(const mat &r);
    rowvec gradientJastrowEvaluate(const mat &r, uint i);
    double alpha,beta;
private:
    double r1, r2, rij,r1r2;
    double eIntEnergy,eContributor;
    double correlation, argument;
    double ddJastrowFactor;
    rowvec dJastrowFactor;
    double JastrowFactor;

protected:
    mat a;

};

#endif // JASTROW_H
