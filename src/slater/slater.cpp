#include "slater.h"

Slater::Slater(const uint &nParticles, Orbitals* orbitals):
    orbitals(orbitals),
    N(nParticles/2),
    DUp(zeros(N,N)),
    DDown(zeros(N,N)),
    DUpNew(zeros(N,N)),
    DDownNew(zeros(N,N)),
    DUpInv(zeros(N,N)),
    DDownInv(zeros(N,N)),
    DUpInvNew(zeros(N,N)),
    DDownInvNew(zeros(N,N)),
    rNew(zeros(nParticles,3))

{
}


/************************************************************
Name:
Description:
*/
void Slater::initializeSD(const mat &r){

    evaluateSD(r);
    DUpNew=DUp;
    DDownNew=DDown;

    DUpInv=inv(DUp);
    DDownInv=inv(DDown);

    DUpInvNew=DUpInv;
    DDownInvNew=DDownInv;

}


/************************************************************
Name:
Description:
*/
double Slater::evaluateSD(const mat &r){
    for (uint p = 0; p < N; p++) {
        for (uint qNum = 0; qNum < N; qNum++) {
            DUp(p,qNum) = orbitals->orbitalEvaluate(r,qNum,p);
            DDown(p, qNum) = orbitals->orbitalEvaluate(r,qNum,p+N);
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
void Slater::updateSlater(){

    if (activeParticle < N) {
        for (uint qNum = 0; qNum  < N; qNum ++) {
            DUpNew(activeParticle, qNum) =orbitals->orbitalEvaluate(rNew,qNum ,activeParticle);
        }
    } else {
        for (uint qNum  = 0; qNum  < N; qNum ++) {
            DDownNew(activeParticle - N, qNum) = orbitals->orbitalEvaluate(rNew,qNum ,activeParticle);
        }
    }
    updateSlaterInverse();
    //            DUpInvNew=inv(DUpNew);
    //            DDownInvNew=inv(DDownNew);
}

/************************************************************
Name:
Description:
*/
void Slater::updateSlaterInverse() {

    R = getSDRatio();

    if (activeParticle < N) {
        uint i = activeParticle;
        S = zeros(1, N);

        for (uint j = 0; j < N; j++){
            for (uint l = 0; l < N; l++){
                S(j) += DUpNew(i, l) * DUpInv(l, j);
            }
        }

        for (uint k = 0; k < N; k++) {
            for (uint j = 0; j < N; j++) {
                if (j != i) {
                    DUpInvNew(k, j) = DUpInv(k, j) - DUpInv(k, i) * S(j) / R;
                } else {
                    DUpInvNew(k, j) = DUpInv(k, i) / R;
                }
            }
        }
    }

    if (activeParticle >= N) {
        uint i = activeParticle - N;
        S = zeros(1, N);

        for (uint j = 0; j < N; j++){
            for (uint l = 0; l < N; l++){
                S(j) += DDownNew(i, l) * DDownInv(l, j);
            }
        }

        for (uint k = 0; k < N; k++) {
            for (uint j = 0; j < N; j++) {
                if (j != i) {
                    DDownInvNew(k, j) = DDownInv(k, j) - DDownInv(k, i) * S(j) / R;
                } else {
                    DDownInvNew(k, j) = DDownInv(k, i) / R;
                }
            }
        }
    }

}


/************************************************************
Name:
Description:
*/
double Slater::getSDRatio(){
    double R = 0;
    uint i = activeParticle;

    if (i < N) {
        for (uint j = 0; j < N; j++) {
            R += DUpNew(i,j) * DUpInv(j, i);
        }
    } else {
        for (uint j = 0; j < N; j++) {
            R += DDownNew(i - N,j) * DDownInv(j, i - N);
        }
    }
    return R;
}



/************************************************************
Name:
Description:
*/
void Slater::acceptMove(){

    if (activeParticle < N) {
        for (uint i = 0; i < N; i++){
            DUp(activeParticle, i) = DUpNew(activeParticle, i);
        }
        DUpInv = DUpInvNew;
    }
    else {
        for (uint i = 0; i < N; i++){
            DDown(activeParticle - N, i) = DDownNew(activeParticle - N, i);
        }
        DDownInv = DDownInvNew;
    }
}


/************************************************************
Name:
Description:
*/
void Slater::rejectMove(){

    if (activeParticle < N) {
        for (uint i = 0; i < N; i++){
            DUpNew(activeParticle, i) = DUp(activeParticle,i);
        }
        DUpInvNew = DUpInv;
    }
    else {
        for (uint i = 0; i < N; i++){
            DDownNew(activeParticle - N, i) = DDown(activeParticle - N, i);
        }
        DDownInvNew = DDownInv;
    }

}



/************************************************************
Name:               gradientSlater
Description:
*/

rowvec Slater::gradientSDEvaluate(const mat &r, uint &i) {
    dSD = zeros(1, r.n_cols);

    if (i < N) {
        for (uint j = 0; j < N; j++) {
            dSD += orbitals->gradientOrbitalEvaluate(r,j,i) * DUpInvNew(j, i);
        }
    } else {
        for (uint j = 0; j < N; j++) {
            dSD += orbitals->gradientOrbitalEvaluate(r,j,i) * DDownInvNew(j, i - N);
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
        for (uint j = 0; j < N; j++) {
            ddSD += orbitals->laplaceOrbitalEvaluate(r,j,i)*DUpInv(j, i);
        }
    } else {
        for (uint j = 0; j < N; j++) {
            ddSD += orbitals->laplaceOrbitalEvaluate(r,j,i)*DDownInv(j, i-N);
        }
    }

    return ddSD;
}



/************************************************************
Name:          laplace
Description:
*/
double Slater::getVariationalDerivate(const mat &r)
{
    //    dVSD = 0;
    //    for (uint i = 0; i < r.n_rows; i++){
    //        if (i < N) {
    //            for (uint j = 0; j < N; j++) {
    //                dVSD += orbitals->getVariationalDerivative(r,j,i)*DUpInv(j, i);
    //            }
    //        } else {
    //            for (uint j = 0; j < N; j++) {
    //                dVSD += orbitals->getVariationalDerivative(r,j,i)*DDownInv(j, i-N);
    //            }
    //        }
    //    }
    //    return dVSD;


    dVSD = 0;
    for (uint i = 0; i < N; i++){
        for (uint j = 0; j < N; j++) {
            dVSD += orbitals->getVariationalDerivative(r,j,i)*DUpInv(j, i)
                    +orbitals->getVariationalDerivative(r,j,i+N)*DDownInv(j, i) ;
        }
    }
    return dVSD;
}

