#ifndef NUMERICALKINETIC_H
#define NUMERICALKINETIC_H

#include "src/Kinetic/kinetic.h"


class NumericalKinetic : public Kinetic
{
public:
    NumericalKinetic();
    double evaluate(int nDimensions, int nParticles,const mat &r);

private:
    mat rPlus ,rMinus;
    double waveFunctionMinus, waveFunctionPlus, waveFunctionCurrent;
    double KineticEnergy,h,h2;
};

#endif // NUMERICALKINETIC_H
