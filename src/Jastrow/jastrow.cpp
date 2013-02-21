#include "jastrow.h"

Jastrow::Jastrow()
{
}

/************************************************************
Name:               JastrowFactor
Description:        computes jastrowfactor
*/
double Jastrow::JastrowExponential(int nParticles, const mat &r){

    correlation=0;
    for (int i=0; i<nParticles-1; i++) {
        for (int j=i+1; j<nParticles; j++) {
            rij= norm( r.row(i)-r.row(j) ,2);

            correlation+=rij/(2+2*beta*rij);
        }
    }

    JastrowFactor = exp(correlation);
    return JastrowFactor;

}


/************************************************************
Name:               LaplaceJastrowEvaluate
Description:
*/

double Jastrow::LaplaceJastrowEvaluate(const mat &r){
    r1 = norm(r.row(0), 2);
    r2 = norm(r.row(1), 2);
    rij = norm(r.row(0) - r.row(1), 2);
    r1r2=dot(r.row(0),r.row(1));


    eIntEnergy= 1./(2 * pow(1+beta*rij,2) );
    eContributor= ((alpha*r1 +alpha*r2)/rij)*(1 - r1r2/(r1*r2) );
    ddJastrowFactor= -2*eIntEnergy*(eContributor-eIntEnergy - 2/rij + 2*beta/(1+beta*rij));
   return ddJastrowFactor;

}
