#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>

using namespace arma;
using namespace std;
using namespace libconfig;


class Potential
{
public:
    Potential();
    virtual double evaluate(const mat &r) = 0;
    void loadConfiguration(Config *cfg);

protected:
 int charge;

};

#endif // POTENTIAL_H
