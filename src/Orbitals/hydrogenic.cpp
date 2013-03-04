#include "hydrogenic.h"

Hydrogenic::Hydrogenic()
{
}


double Hydrogenic::orbitalEvaluate(const mat &r, int qNum, int Particle){
    rNorm= norm(r.row(Particle),2);
    if(qNum==0){
        phi= exp(-k*rNorm);
    }
    else if (qNum==1){
        phi = (k*rNorm - 2)*exp(-0.5*k*rNorm);
    }
    else{
        cerr << "Orbital doesn't exist!"<<endl;
        exit(1);}

    return phi;
}


rowvec Hydrogenic::gradientOrbitalEvaluate(const mat &r, int qNum, int Particle){
    rNorm= norm(r.row(Particle),2);
    dphi=zeros(1,r.n_cols);

    if(qNum==0){
        dphi= (-k/rNorm)*r.row(Particle);//*exp(-k*rNorm);
    }
    else if (qNum==1){
        dphi = (-k/(2*rNorm))*(k*rNorm-4)*r.row(Particle);// *exp(-0.5*k*r.row(Particle));
    }
    else{
        cerr << "Orbital doesn't exist!"<<endl;
        exit(1);
    }

    return dphi;
}


double Hydrogenic::laplaceOrbitalEvaluate(const mat &r, int qNum, int Particle){
    rNorm= norm(r.row(Particle),2);

    if(qNum==0){
        ddphi= k*(k*rNorm-2)/rNorm;
    }
    else if (qNum==1){
        ddphi = k*(k*rNorm-8)*(k*rNorm-2)/(4*rNorm);
    }
    else{
        cerr << "Orbital doesn't exist!"<<endl;
        exit(1);
    }

    return ddphi;
}
