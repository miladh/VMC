#include "hlikewavefunction.h"
#include"src/Jastrow/nojastrow.h"
#include"src/Jastrow/padejastrow.h"
#include"src/Orbitals/hydrogenic.h"


HLikeWavefunction::HLikeWavefunction(const uint &nParticles):
    Wavefunction(nParticles),
    dHydrogenic(zeros(nParticles,3)),
    dJastrow(zeros(nParticles,3))
{
}



/************************************************************
Name:
Description:
*/
void HLikeWavefunction::initializewavefunction(const mat &r)
{
    slater->initializeSD(r);
    jas->initializeJastrow(r);
}

/************************************************************
Name:
Description:
*/
double HLikeWavefunction::wavefunction(const mat &r)
{
    TrialWavefunction =slater->evaluateSD(r);
    TrialWavefunction*=exp(jas->evaluateJastrow(r));

    return TrialWavefunction;
}


/************************************************************
Name:
Description:
*/
void HLikeWavefunction::activeParticle(const mat &r,const uint &i)
{
    slater->setActiveParticleAndCurrentPosition(r,i);
    jas->setActiveParticleAndCurrentPosition(r,i);
}


/************************************************************
Name:
Description:
*/
void HLikeWavefunction::updateWavefunction()
{
    slater->updateSlater();

}


/************************************************************
Name:
Description:
*/
double HLikeWavefunction::getRatio()
{

    return slater->getSDRatio()*jas->getJasRatio();
}


/************************************************************
Name:
Description:
*/
void HLikeWavefunction::rejectMove()
{
    slater->rejectMove();
    jas->rejectMove();
}
/************************************************************
Name:
Description:
*/
void HLikeWavefunction::acceptMove()
{
    slater->acceptMove();
    jas->acceptMove();
}



/************************************************************
Name:          Gradient
Description:
*/
mat HLikeWavefunction::gradient(const mat &r){

    if(useAnalyticGradient){
        for (uint i = 0; i < r.n_rows; i++){
            dHydrogenic.row(i)=slater->gradientSDEvaluate(r,i);
            dJastrow.row(i)=jas->gradientJastrowEvaluate(r,i);
        }
        dwavefunction=dHydrogenic+dJastrow;

        return dwavefunction;
    }
    else{
        return gradientNumerical(r);
    }

}


/************************************************************
Name:          laplace
Description:
*/
double HLikeWavefunction::laplace(const mat &r){

    if(useAnalyticLaplace){
        ddwavefunction = 0;
        for (uint i = 0; i < r.n_rows; i++){

            ddwavefunction+=slater->laplaceSDEvaluate(r,i)+
                    2*dot(slater->gradientSDEvaluate(r,i),jas->gradientJastrowEvaluate(r,i));
        }

        ddwavefunction+= jas->laplaceJastrowEvaluate(r);

        return ddwavefunction;
    }
    else{

        return laplaceNumerical(r);
    }

}

/************************************************************
Name:          laplace
Description:
*/
vec HLikeWavefunction::getVariationalDerivate(const mat &r)
{
    dvariational(0) = slater->getVariationalDerivate(r);
    dvariational(1) = jas->getVariationalDerivative(r);
    return dvariational;
}


















