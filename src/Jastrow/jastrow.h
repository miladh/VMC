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

protected:
    double ddJastrowFactor;
    double JastrowFactor;
    double deltaJastrow;
    uint activeParticle;
    rowvec dJastrowFactor;

public:
    Jastrow();
    virtual double evaluateJastrow(const mat& r)=0;
    virtual double laplaceJastrowEvaluate(const mat& r)=0;
    virtual rowvec gradientJastrowEvaluate(const mat& r, uint i)=0;
    virtual double getJasRatio()=0;
    virtual void initializeJastrow(const mat& r)=0;
    virtual void setActiveParticleAndCurrentPosition(const mat& r, const uint& i)=0;
    virtual void acceptMove()=0;
    virtual void rejectMove()=0;
    virtual double getVariationalDerivative(const mat &r) = 0;
};

#endif // JASTROW_H
