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
        phi = (k*rNorm-2 )*exp(-0.5*k*rNorm);
    }
    else if (qNum == 2) {
        phi = r(Particle,2)*exp(-0.5*k*rNorm);
    }
    else if (qNum == 3) {
        phi =r(Particle,0)*exp(-0.5*k*rNorm);
    }
    else if (qNum == 4) {
        phi =  r(Particle,1)*exp(-0.5*k*rNorm);
    }
    else if(qNum==5){
        phi= (1 - 2*k*rNorm/3 + 2*k*k*rNorm*rNorm/27)*
             exp(-k*rNorm/3);
    }
    else{
        cerr << "Orbital doesn't exist!"<<endl;
        exit(1);}

    return phi;
}


rowvec Hydrogenic::gradientOrbitalEvaluate(const mat &r, int qNum, int Particle){
    rNorm= norm(r.row(Particle),2);
    dphi=zeros(1,r.n_cols);

    // 1s orbital
    if(qNum==0){
        dphi= (-k/rNorm)*r.row(Particle)*exp(-k*rNorm);
    }

    // 2s orbital
    else if (qNum==1){
        dphi = (-k/(2*rNorm))*(k*rNorm-4)*r.row(Particle)*exp(-0.5*k*rNorm);
    }

    // 2p orbital
    else if (qNum == 2) {
        dphi(0,0) = r(Particle,0)*r(Particle,2);
        dphi(0,1) = r(Particle,1)*r(Particle,2);
        dphi(0,2) = r(Particle,2)*r(Particle,2)-2*rNorm/k;
        dphi*= (-k/(2*rNorm))*exp(-0.5*k*rNorm);
    }
    else if (qNum == 3) {
        dphi(0,0) = r(Particle,0)*r(Particle,0)-2*rNorm/k;
        dphi(0,1) = r(Particle,0)*r(Particle,1);
        dphi(0,2) = r(Particle,0)*r(Particle,2);
        dphi*= (-k/(2*rNorm))*exp(-0.5*k*rNorm);
    }
    else if (qNum == 4) {
        dphi(0,0) = r(Particle,1)*r(Particle,0);
        dphi(0,1) = r(Particle,1)*r(Particle,1)-2*rNorm/k;
        dphi(0,2) = r(Particle,1)*r(Particle,2);
        dphi*= (-k/(2*rNorm))*exp(-0.5*k*rNorm);
    }

    // 3s orbital
    else if(qNum==5){
        dphi=r.row(Particle)*(-1 + 10*k*rNorm/27 - 2*k*k*rNorm*rNorm/81)
                *k*exp(-k*rNorm/3);
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
        ddphi= k*(k*rNorm-2)/rNorm*exp(-k*rNorm);
    }
    else if (qNum==1){
        ddphi = k*(k*rNorm-8)*(k*rNorm-2)/(4*rNorm)*exp(-0.5*k*rNorm);
    }
    else if (qNum==2){
        ddphi = k*(k*rNorm-8)/(4*rNorm)*
                r(Particle,2)*exp(-0.5*k*rNorm);
    }
    else if (qNum==3){
        ddphi = k*(k*rNorm-8)/(4*rNorm)*
                r(Particle,0)*exp(-0.5*k*rNorm);
    }
    else if (qNum==4){
        ddphi = k*(k*rNorm-8)/(4*rNorm)*
                r(Particle,1)*exp(-0.5*k*rNorm);
    }
    else if(qNum==5){
        ddphi= (-2 + 13*k*rNorm/9 - 2*k*k*rNorm*rNorm/9
               +2*k*k*k*rNorm*rNorm*rNorm/243)*k/rNorm
              *exp(-k*rNorm/3);
    }
    else{
        cerr << "Orbital doesn't exist!"<<endl;
        exit(1);
    }

    return ddphi;
}
