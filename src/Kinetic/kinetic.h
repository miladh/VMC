#ifndef KINETIC_H
#define KINETIC_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include "src/Wavefunction/wavefunction.h"

using namespace arma;
using namespace std;
using namespace libconfig;


class Kinetic
{
public:
    Kinetic();
     virtual double evaluate(int nDimensions, int nParticles,const mat &r) = 0;

    Wavefunction* wf;
    double alpha,beta;
};

#endif // KINETIC_H