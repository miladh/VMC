#ifndef NOJASTROW_H
#define NOJASTROW_H
#include"jastrow.h"

class NoJastrow : public Jastrow
{
public:
    NoJastrow();

    inline double evaluateJastrow(const mat&){return 0;}
    inline rowvec gradientJastrowEvaluate(const mat&, uint){ return zeros(1,3);}
    inline double laplaceJastrowEvaluate(const mat&){return 0;}
    inline double getJasRatio(){return 1;}
    inline void initializeJastrow(const mat&){}
    inline void setActiveParticleAndCurrentPosition(const mat&, const uint&){}
    inline void acceptMove(){}
    inline void rejectMove(){}
    inline double getVariationalDerivative(const mat&){return 0;}

};

#endif // NOJASTROW_H
