#include "wavefunction.h"


Wavefunction::Wavefunction(Config* cfg, Orbitals* orbitals, Jastrow* jastrow):
    cfg(cfg),
    orbitals(orbitals),
    jas(jastrow)
{
        loadAndSetConfiguration();
}

//************************************************************
void Wavefunction::initializeWavefunction(const mat &r)
{
    slater->initializeSD(r);
    jas->initializeJastrow(r);
}

//************************************************************
double Wavefunction::evaluateWavefunction(const mat &r)
{
    trialWavefunction  = slater->evaluateSD(r);
    trialWavefunction *= exp(jas->evaluateJastrow(r));

    return trialWavefunction;
}


//************************************************************
void Wavefunction::activeParticle(const mat &r,const uint &i)
{
    slater->setActiveParticleAndCurrentPosition(r,i);
    jas->setActiveParticleAndCurrentPosition(r,i);
}

//************************************************************
void Wavefunction::updateWavefunction()
{
    slater->updateSlater();

}

//************************************************************
double Wavefunction::getRatio()
{

    return slater->getSDRatio()*jas->getJasRatio();
}

//************************************************************
void Wavefunction::rejectMove()
{
    slater->rejectMove();
    jas->rejectMove();
}

//************************************************************
void Wavefunction::acceptMove()
{
    slater->acceptMove();
    jas->acceptMove();
}



//************************************************************
mat Wavefunction::gradient(const mat &r){

    if(useAnalyticGradient){
        for (uint i = 0; i < r.n_rows; i++){
            dSlater.row(i)=slater->gradientSDEvaluate(r,i);
            dJastrow.row(i)=jas->gradientJastrowEvaluate(r,i);
        }
        dwavefunction = dSlater+dJastrow;

        return dwavefunction;
    }
    else{
        return gradientNumerical(r);
    }

}

//************************************************************
double Wavefunction::laplace(const mat &r){

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

//************************************************************
vec Wavefunction::getVariationalDerivate(const mat &r)
{
    dvariational(0) = slater->getVariationalDerivate(r);
    dvariational(1) = jas->getVariationalDerivative(r);
    return dvariational;
}


//************************************************************
double Wavefunction::laplaceNumerical(const mat &r){

    rPlus = rMinus = r;
    wavefunctionCurrent = evaluateWavefunction(r);
    ddwavefunction = 0;

    for(uint i = 0; i <r.n_rows; i++) {
        for(uint j=0; j <r.n_cols ; j++){
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            wavefunctionMinus = evaluateWavefunction(rMinus);
            wavefunctionPlus = evaluateWavefunction(rPlus);
            ddwavefunction += (wavefunctionMinus + wavefunctionPlus - 2 * wavefunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j)= r(i,j);
        }
    }

    ddwavefunction =h2 * ddwavefunction / wavefunctionCurrent;

    return ddwavefunction;
}



//************************************************************
mat Wavefunction::gradientNumerical(const mat &r){

    rPlus = rMinus = r;

    wavefunctionCurrent = evaluateWavefunction(r);
    for(uint i =0; i<r.n_rows ; i++){
        for(uint j =0; j<r.n_cols ; j++){
            rPlus(i,j)+=hGrad;
            rMinus(i,j)-=hGrad;
            wavefunctionMinus = evaluateWavefunction(rMinus);
            wavefunctionPlus = evaluateWavefunction(rPlus);
            dwavefunction(i,j)=(wavefunctionPlus-wavefunctionMinus)/(2*wavefunctionCurrent*hGrad);
            rPlus(i,j)=r(i,j);
            rMinus(i,j)=r(i,j);
        }
    }
    return dwavefunction;

}


//************************************************************
void Wavefunction::loadAndSetConfiguration(){
    nParticles  = cfg->lookup("setup.nParticles");
    nDimensions = cfg->lookup("setup.nDimensions");
    useAnalyticLaplace  = cfg->lookup("AppSettings.useAnalyticLaplace");
    useAnalyticGradient = cfg->lookup("AppSettings.useAnalyticGradient");

    dSlater       = zeros(nParticles,nDimensions);
    dJastrow      = zeros(nParticles,nDimensions);
    dwavefunction = zeros<mat>(nParticles,nDimensions);
    rPlus         = zeros<mat>(nParticles, nDimensions);
    rMinus        = zeros<mat>(nParticles, nDimensions);
    dvariational  = zeros<vec>(2);

    slater  = new Slater(nParticles,orbitals);

}





