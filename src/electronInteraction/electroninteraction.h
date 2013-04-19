#ifndef ELECTRONINTERACTION_H
#define ELECTRONINTERACTION_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>

using namespace arma;
using namespace std;
using namespace libconfig;

class ElectronInteraction
{
public:
    ElectronInteraction();

    virtual double evaluate(const mat& r)=0;
};

#endif // ELECTRONINTERACTION_H
