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
    Kinetic(Config *cfg);
    double evaluate(int nParticles,const mat &r);
    double KineticEnergy;
    Wavefunction* wf;

private:
    double ddwaveFunction;
    Config *cfg;

};

#endif // KINETIC_H
