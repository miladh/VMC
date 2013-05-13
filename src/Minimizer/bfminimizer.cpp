#include <src/Minimizer/bfminimizer.h>


BFMinimizer::BFMinimizer(Config *cfg, const int& myRank, const int& nProcess):
    Minimizer(cfg,myRank,nProcess)
{
    loadAndSetConfiguration();
}

//*****************************************************************************
void BFMinimizer::runMinimizer(){

    myfile.open ("../vmc/DATA/results");
    myfile << "Alpha   " << "Beta    " << "Energy       " << "Variance    "
           << "Sigma       "<< "Acceptance    "<< endl;

    for(int i = 0; i < nVarAlpha; i++){
        parser->alpha = alpha;
        beta          = minBeta;

        for(int j = 0; j < nVarBeta; j++){
            parser->beta = beta;
            parser->setup();
            getResultsAndWrite();
            beta += stepBeta;
        }
        alpha += stepAlpha;
    }

    myfile.close();
}

//*****************************************************************************
void BFMinimizer::getResultsAndWrite()
{

    energy = parser->getEnergy();
    energySquared = parser->getEnergySquared();
    variance = parser->getVariance();
    sigma = parser->getSigma();
    acceptance= parser->getAcceptanceRate();

    if (myRank == 0) {
    myfile << alpha <<"     "<<  beta <<"     "<<energy
           <<"     "<< variance <<"     "<<sigma
           <<"     "<< acceptance<< endl;
    }
}

//*****************************************************************************
void BFMinimizer::loadAndSetConfiguration()
{
    minAlpha  = cfg->lookup("setup.MinimizerSettings.BFMinSettings.minalpha");
    maxAlpha  = cfg->lookup("setup.MinimizerSettings.BFMinSettings.maxalpha");
    minBeta   = cfg->lookup("setup.MinimizerSettings.BFMinSettings.minbeta");
    maxBeta   = cfg->lookup("setup.MinimizerSettings.BFMinSettings.maxbeta");
    nVarAlpha = cfg->lookup("setup.MinimizerSettings.BFMinSettings.nVarAlpha");
    nVarBeta  = cfg->lookup("setup.MinimizerSettings.BFMinSettings.nVarBeta");

    alpha     = minAlpha;
    beta      = minBeta;
    stepAlpha = (maxAlpha-minAlpha)/nVarAlpha;
    stepBeta  = (maxBeta-minBeta)/nVarBeta;

}


