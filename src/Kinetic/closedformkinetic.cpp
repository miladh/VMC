#include "closedformkinetic.h"

ClosedFormKinetic::ClosedFormKinetic():
    charge(2)
{
}

double ClosedFormKinetic::evaluate(int nDimensions, int nParticles,const mat &r){


    r1 = norm(r.row(0), 2);
    r2 = norm(r.row(1), 2);
    rij = norm(r.row(0) - r.row(1), 2);
    r1r2=dot(r.row(0),r.row(1));


    E_L1 = (alpha - charge)*(1.0/r1 + 1.0/r2) + 1.0/rij - alpha*alpha;
    eIntEnergy= 1./(2 * pow(1+beta*rij,2) );
    eContributor= ((alpha*r1 +alpha*r2)/rij)*(1 - r1r2/(r1*r2) );
    localEnergy= E_L1+ eIntEnergy*(eContributor-eIntEnergy - 2/rij + 2*beta/(1+beta*rij));

    return localEnergy;
}


