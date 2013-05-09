#ifndef PADEJASTROW_H
#define PADEJASTROW_H


#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include "src/Jastrow/jastrow.h"

using namespace arma;
using namespace std;
using namespace libconfig;

class PadeJastrow :public Jastrow
{
public:
    PadeJastrow(const uint nParticles, const double &beta);
    double evaluateJastrow(const mat &r);
    double laplaceJastrowEvaluate(const mat &r);
    rowvec gradientJastrowEvaluate(const mat &r, uint i);
    double getJasRatio();
    void initializeJastrow(const mat &r);
    void setActiveParticleAndCurrentPosition(const mat &r, const uint &i);
    void acceptMove();
    void rejectMove();
    double getVariationalDerivative(const mat &r);


private:
    uint nParticles;
    double beta;
    double rij;
    double correlation;
    double r_ki,b_ij;
    mat a, rOld,rNew;

    void setaValues();


};

#endif // PADEJASTROW_H
