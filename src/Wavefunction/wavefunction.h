#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include <src/slater/slater.h>
#include <src/Jastrow/jastrow.h>
#include <src/includes/Defines.h>

using namespace arma;
using namespace std;
using namespace libconfig;

class Wavefunction
{

public:
    Wavefunction(Config *cfg, Orbitals *orbitals, Jastrow *jastrow);


    void activeParticle(const mat &r,const uint &i);
    double evaluateWavefunction(const mat &r);
    void initializeWavefunction(const mat &r);
    void updateWavefunction();
    void acceptMove();
    void rejectMove();
    double getRatio();
    double laplace(const mat &r);
    mat gradient(const mat &r);
    vec getVariationalDerivate(const mat &r);


private:
    uint nParticles, nDimensions;
    double trialWavefunction;
    double ddwavefunction;
    double wavefunctionMinus, wavefunctionPlus, wavefunctionCurrent;
    bool useAnalyticGradient,useAnalyticLaplace;
    mat dwavefunction,dSlater,dJastrow;
    mat rPlus ,rMinus;
    vec dvariational; 

    Config* cfg;
    Orbitals* orbitals;
    Jastrow* jas;
    Slater *slater;

    void loadAndSetConfiguration();
    double laplaceNumerical(const mat &r);
    mat gradientNumerical(const mat &r);


};

#endif // WAVEFUNCTION_H
