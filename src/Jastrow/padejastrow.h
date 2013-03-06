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
    PadeJastrow(const uint nParticles);
    void setaValues(const uint &nParticles);
    double evaluateJastrow(const mat &r);
    double laplaceJastrowEvaluate(const mat &r);
    rowvec gradientJastrowEvaluate(const mat &r, uint i);

private:
    double rij;
    double correlation;
    double r_ki,b_ij;


};

#endif // PADEJASTROW_H
