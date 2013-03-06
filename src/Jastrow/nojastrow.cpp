#include "nojastrow.h"

NoJastrow::NoJastrow(const uint nParticles):
    Jastrow(nParticles)
{
}


void NoJastrow::setaValues(const uint &nParticles){
}


double NoJastrow::evaluateJastrow(const mat &r){
    return 1;
}


rowvec NoJastrow::gradientJastrowEvaluate(const mat &r, uint i) {
    return zeros(1,r.n_cols);
}


double NoJastrow::laplaceJastrowEvaluate(const mat &r){
    return 0;
}
