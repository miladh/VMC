#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include <src/VMCApp/vmcapp.h>

using namespace arma;
using namespace std;
using namespace libconfig;


class Minimizer
{
public:
    Minimizer(const int &, const int &);

    virtual void runMinimizer()=0;
    virtual void loadConfiguration(Config *cfg)=0;
};

#endif // MINIMIZER_H
