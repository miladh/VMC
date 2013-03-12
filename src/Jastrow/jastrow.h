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
    uint nParticles;
    mat a;
    double ddJastrowFactor;
    rowvec dJastrowFactor;
    double JastrowFactor;
    uint activeParticle;
    double deltaJastrow;

public:
    Jastrow(const uint &nParticles);
    virtual void setaValues(const uint &nParticles)=0;
    virtual double evaluateJastrow(const mat &r)=0;
    virtual double laplaceJastrowEvaluate(const mat &r)=0;
    virtual rowvec gradientJastrowEvaluate(const mat &r, uint i)=0;
    virtual double getJasRatio()=0;
    virtual void initializeJastrow(const mat &r)=0;
    virtual void setActiveParticleAndCurrentPosition(const mat &r, const uint &i)=0;
    virtual void acceptMove()=0;
    virtual void rejectMove()=0;

    double alpha,beta;
    mat rOld,rNew;
};

#endif // JASTROW_H
