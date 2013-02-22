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
    double evaluate(const mat &r);
    Wavefunction* wf;

private:
    double KineticEnergy;
    double ddwavefunction;

};

#endif // KINETIC_H
