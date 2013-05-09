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
    Kinetic(Wavefunction *wavefunction);
    double evaluate(const mat &r);


private:
    Wavefunction* wavefunction;
};

#endif // KINETIC_H
