#include <src/Minimizer/steepestdescent.h>

SteepestDescent::SteepestDescent(Config *cfg, const int &myRank, const int &nProcess):
    Minimizer(cfg,myRank,nProcess)
{
    loadAndSetConfiguration();
}

//*****************************************************************************
void SteepestDescent::runMinimizer(){

    vec old = zeros<vec>(2);

    myfile.open ("../vmc/DATA/results");
    myfile << "Alpha   " << "Beta    " << "Energy       " << "Variance    "
           << "Sigma       "<< "Acceptance    "<< endl;


    for(int i=0; i < 100; i++){

        parser->alpha = alpha;
        parser->beta  = beta;
        parser->setup();

        variationalDerivate = parser->getVariationalDerivate();
        getResultsAndWrite();

        if(old(0)/variationalDerivate(0) >= 0){
            step *=1.25;
            alpha -= step*signFunc(variationalDerivate(0));
        }else{
            step*=0.5;
            alpha -= step*signFunc(variationalDerivate(0));
        }


        if(old(1)/variationalDerivate(1) >= 0){
            step*=1.25;
            beta -= step*signFunc(variationalDerivate(1));
        }else{
            step*=0.5;
            beta -= step*signFunc(variationalDerivate(1));
        }

        old = variationalDerivate;
    }
        myfile.close();
}

//*****************************************************************************
int SteepestDescent::signFunc(double varDer){

    if(varDer < 0){
        return -1;
    }else{
        return 1;
    }
}


//*****************************************************************************
void SteepestDescent::getResultsAndWrite()
{
    energy = parser->getEnergy();
    energySquared = parser->getEnergySquared();
    variance=parser->getVariance();
    sigma=parser->getSigma();
    acceptance= parser->getAcceptanceRate();

    if (myRank == 0) {
    myfile << alpha <<"     "<<  beta <<"     "<<energy
           <<"     "<< variance <<"     "<<sigma
           <<"     "<< acceptance<< endl;
    }
}

//*****************************************************************************
void SteepestDescent::loadAndSetConfiguration()
{
    alpha = cfg->lookup("setup.MinimizerSettings.SDMinSettings.InitAlpha");
    beta  = cfg->lookup("setup.MinimizerSettings.SDMinSettings.InitBeta");
    step  = cfg->lookup("setup.MinimizerSettings.SDMinSettings.step");
    variationalDerivate = zeros<vec>(2);
}
