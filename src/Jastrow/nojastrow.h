#ifndef NOJASTROW_H
#define NOJASTROW_H
#include"jastrow.h"

class NoJastrow : public Jastrow
{
public:
    NoJastrow(const uint nParticles);

    void setaValues(const uint &nParticles);
    double evaluateJastrow(const mat &r);
    double laplaceJastrowEvaluate(const mat &r);
    rowvec gradientJastrowEvaluate(const mat &r, uint i);
};

#endif // NOJASTROW_H
