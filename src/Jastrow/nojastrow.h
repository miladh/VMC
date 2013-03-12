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
    double getJasRatio();
    void initializeJastrow(const mat &r);
    void setActiveParticleAndCurrentPosition(const mat &r, const uint &i);
    void acceptMove();
    void rejectMove();
};

#endif // NOJASTROW_H
