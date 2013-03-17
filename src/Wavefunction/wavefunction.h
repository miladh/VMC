#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include "src/Orbitals/hydrogenic.h"
#include "src/slater/slater.h"
#include "src/Jastrow/jastrow.h"

using namespace arma;
using namespace std;
using namespace libconfig;

class Wavefunction
{
protected:
    uint nParticles;
    bool useAnalyticGradient,useAnalyticLaplace;
    double TrialWavefunction;
    double kineticEnergy;
    double ddwavefunction;
    mat dwavefunction;



private:
    double hGrad,h,h2;
    mat rPlus ,rMinus;
    double wavefunctionMinus, wavefunctionPlus, wavefunctionCurrent;

public:
    Wavefunction(const uint &nParticles);

    virtual double wavefunction(const mat &r) = 0;
    virtual double laplaceNumerical(const mat &r);
    virtual mat gradientNumerical(const mat &r);

    virtual double laplace(const mat &r){
        return laplaceNumerical(r);}

    virtual mat gradient(const mat &r){
        return gradientNumerical(r);}

    void loadConfiguration(Config *cfg);

    Slater *slater;
    Jastrow* jas;


    virtual void activeParticle(const mat &r,const uint &i)=0;
    virtual void updateWavefunction()=0;
    virtual double getRatio()=0;
    virtual void acceptMove()=0;
    virtual void rejectMove()=0;
    virtual void initializewavefunction(const mat &r)=0;


};

#endif // WAVEFUNCTION_H
