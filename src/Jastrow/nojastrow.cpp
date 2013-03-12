#include "nojastrow.h"

NoJastrow::NoJastrow(const uint nParticles):
    Jastrow(nParticles)
{
}


void NoJastrow::setaValues(const uint &nParticles){
}


double NoJastrow::evaluateJastrow(const mat &r){
    return 0;
}


rowvec NoJastrow::gradientJastrowEvaluate(const mat &r, uint i) {
    return zeros(1,r.n_cols);
}


double NoJastrow::laplaceJastrowEvaluate(const mat &r){
    return 0;
}


double NoJastrow::getJasRatio(){
    return 1;
}

void NoJastrow::initializeJastrow(const mat &r){
}


void NoJastrow::setActiveParticleAndCurrentPosition(const mat &r, const uint &i){
}


void NoJastrow::acceptMove(){
}

void NoJastrow::rejectMove(){
}
