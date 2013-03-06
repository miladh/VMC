#include "slater.h"
#include "src/Orbitals/hydrogenic.h"

Slater::Slater(const uint &nParticles):
    N(nParticles/2),
    orbitals(new Hydrogenic),
    DUp(zeros(N,N)),
    DDown(zeros(N,N)),
    DUpNew(zeros(N,N)),
    DDownNew(zeros(N,N)),
    DUpInv(zeros(N,N)),
    DDownInv(zeros(N,N)),
    DUpInvNew(zeros(N,N)),
    DDownInvNew(zeros(N,N)),
    rNew(zeros(nParticles,nParticles))

{
}


/************************************************************
Name:
Description:
*/
double Slater::initializeSD(const mat &r){

    initialSD=evaluateSD(r);

    DUpNew=DUp;
    DDownNew=DDown;

    DUpInv=inv(DUp).t();
    DDownInv=inv(DDown).t();
    DUpInvNew=DUpInv;
    DDownInvNew=DDownInv;

    return initialSD;

}


/************************************************************
Name:
Description:
*/
double Slater::evaluateSD(const mat &r){
    for (uint i = 0; i < N; i++) {
        for (uint qNum = 0; qNum < N; qNum++) {
            DUp(qNum, i) = orbitals->orbitalEvaluate(r,qNum,i);
            DDown(qNum, i) = orbitals->orbitalEvaluate(r,qNum,i+N);
        }
    }

    return det(DUp)*det(DDown);

}

/************************************************************
Name:
Description:
*/
void Slater::setActiveParticleAndCurrentPosition(const mat &r, uint i ){
    activeParticle=i;
    rNew=r;
}

/************************************************************
Name:
Description:
*/
double Slater::updateSlater(){

    if (activeParticle < N) {
         for (uint qNum = 0; qNum  < N; qNum ++) {
             DUpNew(qNum , activeParticle) =orbitals->orbitalEvaluate(rNew,qNum ,activeParticle);
         }
     } else {
         for (uint qNum  = 0; qNum  < N; qNum ++) {
             DDownNew(qNum , activeParticle - N) = orbitals->orbitalEvaluate(rNew,qNum ,activeParticle);
         }
     }

    DUpInvNew=inv(DUpNew).t();
    DDownInvNew=inv(DDownNew).t();

    return det(DUpNew)*det(DDownNew);
}


/************************************************************
Name:
Description:
*/
double Slater::getSDRatio(){
    double R = 0;
    uint i = activeParticle;

     if (i < N) { // Spin up
         for (uint j = 0; j < N; j++) {
             R += DUpNew(j, i) * DUpInv(j, i);
         }
     } else {
         for (uint j = 0; j < N; j++) {
             R += DDownNew(j, i - N) * DDownInv(j, i - N);
         }
     }
     return R;
}



/************************************************************
Name:
Description:
*/
void Slater::acceptNewPosition(){
    if (activeParticle < N) {
          for (uint i = 0; i < N; i++)
              DUp(i, activeParticle) = DUpNew(i, activeParticle);
          DUpInv = DUpInvNew;
      } else {
          for (uint i = 0; i < N; i++)
              DDown(i, activeParticle - N) = DDownNew(i, activeParticle - N);
          DDownInv = DDownInvNew;
      }
}






/************************************************************
Name:               gradientSlater
Description:
*/

rowvec Slater::gradientSDEvaluate(const mat &r, uint &i) {
    dSD= zeros(1, r.n_cols);

    if (i < N) {
        for (uint qNum = 0; qNum < N; qNum++) { // Spin up.
            dSD += orbitals->gradientOrbitalEvaluate(r,qNum,i)* DUpInv(qNum, i);
        }
    } else {
        for (uint qNum = 0; qNum < N; qNum++) { // Spin down.
            dSD += orbitals->gradientOrbitalEvaluate(r,qNum,i) * DDownInv(qNum, i - N);
        }
    }
    return dSD;
}


/************************************************************
Name:               laplaceSlater
Description:
*/

double Slater::laplaceSDEvaluate(const mat &r,const uint &i) {

    ddSD = 0;
    if (i < N) {
        for (uint qNum = 0; qNum < N; qNum++) { // Spin up.
            ddSD += orbitals->laplaceOrbitalEvaluate(r,qNum,i)*DUpInv(qNum, i);
        }
    } else {
        for (uint qNum = 0; qNum < N; qNum++) { // Spin down.
            ddSD += orbitals->laplaceOrbitalEvaluate(r,qNum,i)*DUpInv(qNum, i-N);
        }
    }

    return ddSD;
}


