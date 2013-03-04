#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include"src/Orbitals/hydrogenic.h"
#include "src/slater/slater.h"
#include "src/Jastrow/jastrow.h"

using namespace arma;
using namespace std;
using namespace libconfig;

class Wavefunction
{
public:
    Wavefunction();

    virtual double wavefunction(const mat &r) = 0;
    virtual double laplaceNumerical(const mat &r);
    virtual mat gradientNumerical(const mat &r);

    virtual double laplace(const mat &r){
        return laplaceNumerical(r);}

    virtual mat gradient(const mat &r){
        return gradientNumerical(r);}

    void loadConfiguration(Config *cfg);

    Slater slaterDet;
    Orbitals* orbitals;
    Jastrow jas;



protected:
    bool useAnalyticGradient,useAnalyticLaplace;
    double TrialWavefunction;
    double kineticEnergy;
    double ddwavefunction;
    mat dwavefunction;



private:
    mat rPlus ,rMinus;
    rowvec hVec;
    double wavefunctionMinus, wavefunctionPlus, wavefunctionCurrent;
    double hGrad,h,h2;

};

#endif // WAVEFUNCTION_H
