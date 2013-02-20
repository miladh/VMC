#include "jastrowwavefunction.h"

JastrowWavefunction::JastrowWavefunction()
{
}

/************************************************************
Name:               JastrowWavefunction
Description:        jastrow wavefunction wavefunction
*/

double JastrowWavefunction::waveFunction(int nParticles, const mat &r)
{
    argument=0.0;

    for (int i=0; i<nParticles; i++) {
        rSingleParticle = norm(r.row(i),2);
        argument += rSingleParticle;
    }

    TrialWaveFunction = exp(-alpha* argument);

    TrialWaveFunction *=jastrowFactor(nParticles,r);
    return TrialWaveFunction;

}


/************************************************************
Name:               JastrowFactor
Description:        computes jastrowfactor
*/
double JastrowWavefunction::jastrowFactor(int nParticles, const mat &r){
    correlation=0;

    for (int i=0; i<nParticles-1; i++) {
        for (int j=i+1; j<nParticles; j++) {
            rij= norm( r.row(i)-r.row(j) ,2);

            correlation+=rij/(2+2*beta*rij);
        }
    }

    double Jastrow = exp(correlation);
    return Jastrow;

}

/************************************************************
Name:          laplace
Description:
*/
double JastrowWavefunction::laplace(int nParticles, const mat &r, Config* cfg){

    analytic= cfg->lookup("AppSettings.useAnalyticLaplace");

    if(analytic){
    r1 = norm(r.row(0), 2);
    r2 = norm(r.row(1), 2);
    rij = norm(r.row(0) - r.row(1), 2);
    r1r2=dot(r.row(0),r.row(1));

    E_L1 = alpha*(1.0/r1 + 1.0/r2) - alpha*alpha;
    eIntEnergy= 1./(2 * pow(1+beta*rij,2) );
    eContributor= ((alpha*r1 +alpha*r2)/rij)*(1 - r1r2/(r1*r2) );
    ddwaveFunction= 2*alpha*(alpha-1/r1-1/r2)-2*eIntEnergy*(eContributor-eIntEnergy - 2/rij + 2*beta/(1+beta*rij));

   return ddwaveFunction;

    }else{
    return laplaceNumerical(nParticles,r,cfg);
}

}





