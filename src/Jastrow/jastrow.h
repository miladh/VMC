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
    Jastrow(const uint nParticles);
    virtual void setaValues(const uint &nParticles)=0;
    virtual double evaluateJastrow(const mat &r)=0;
    virtual double laplaceJastrowEvaluate(const mat &r)=0;
    virtual rowvec gradientJastrowEvaluate(const mat &r, uint i)=0;
    double alpha,beta;

protected:
    uint nParticles;
    mat a;
    double ddJastrowFactor;
    rowvec dJastrowFactor;
    double JastrowFactor;
};

#endif // JASTROW_H
