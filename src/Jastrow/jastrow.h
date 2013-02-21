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
    double JastrowExponential(int nParticles, const mat &r);
    double LaplaceJastrowEvaluate(const mat &r);
    double alpha,beta;
private:
    double r1, r2, rij,r1r2;
    double eIntEnergy,eContributor;
    double correlation, argument;
    double ddJastrowFactor;
    double JastrowFactor;


};

#endif // JASTROW_H
