#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include"src/Orbitals/hydrogenic.h"
#include "src/Jastrow/jastrow.h"

using namespace arma;
using namespace std;
using namespace libconfig;

class Wavefunction
{
public:
    Wavefunction();
    virtual double waveFunction(int nParticles,const mat &r) = 0;

    virtual double laplace(int nParticles, const mat &r,Config* cfg){
        return laplaceNumerical(nParticles,r,cfg);
    }
//    virtual double gradient(int nParticles, const mat &r){
//        return gradienteNumerical;
//    }


    virtual double laplaceNumerical(int nParticles,const mat &r,Config* cfg);
    //virtual double gradientNumerical(int nParticles,const mat &r);

    double ddwaveFunction;
    double TrialWaveFunction;
    double kineticEnergy;
    bool analytic;
    Config* cfg;
    Orbitals* orbitals;
    Jastrow jas;

private:
    mat rPlus ,rMinus;
    rowvec hVec;
    double waveFunctionMinus, waveFunctionPlus, waveFunctionCurrent;
    double h,h2;



};

#endif // WAVEFUNCTION_H
