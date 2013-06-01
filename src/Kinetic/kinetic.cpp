#include "kinetic.h"

Kinetic::Kinetic(Wavefunction* wavefunction):
    wavefunction(wavefunction)
{
}
//************************************************************
double Kinetic::evaluate(const mat &r){

    return -0.5 * wavefunction->laplace(r);
}
