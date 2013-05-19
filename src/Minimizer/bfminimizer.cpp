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


    if(diatomicSystem){
        Setting& RVector = cfg->lookup("setup.MinimizerSettings.BFMinSettings.R");
        vec RInterval = zeros(3);

        for(int i =0; i < 3; i++){
            RInterval[i] = RVector[i];
        }

        vec RValues = linspace<vec>(RInterval[0],RInterval[1],RInterval[2]);
        parameters.resize(3);


        for(uint i = 0; i < RValues.n_elem; i++){
            parameters[2] = RValues[i];
            minimize();
        }

    }else{
        minimize();
    }

    myfile.close();
}


//*****************************************************************************
void BFMinimizer::minimize()
{

    for(uint i = 0; i < alphaValues.n_elem; i++){
        parameters[0] = alphaValues[i];

        for(uint j = 0; j < betaValues.n_elem; j++){
            parameters[1] = betaValues[j];

            parser->setVariationalParameters(parameters);
            parser->runSolver();
            getResultsAndWrite();
        }
    }

}



//*****************************************************************************
void BFMinimizer::getResultsAndWrite()
{

    energy        = parser->getEnergy();
    energySquared = parser->getEnergySquared();
    variance      = parser->getVariance();
    sigma         = parser->getSigma();
    acceptance    = parser->getAcceptanceRate();

    if (myRank == 0) {
        myfile <<  parameters[0] <<"     " << parameters[1]  << "     "
               <<  parameters[2] <<"     "  << energy <<"     "<< variance
                                 <<"     " << sigma   <<"     "<< acceptance
                                 << endl;
    }
}

//*****************************************************************************
void BFMinimizer::loadAndSetConfiguration()
{
    diatomicSystem = cfg->lookup("setup.system");
    Setting& alphaVector = cfg->lookup("setup.MinimizerSettings.BFMinSettings.alpha");
    Setting& betaVector = cfg->lookup("setup.MinimizerSettings.BFMinSettings.beta");

    vec alphaInterval = zeros(3);
    vec betaInterval = zeros(3);

    for(uint i =0; i < 3; i++){
        alphaInterval[i] = alphaVector[i];
        betaInterval[i]  = betaVector[i];
    }

    parameters.resize(2);
    alphaValues = linspace<vec>(alphaInterval[0],alphaInterval[1],alphaInterval[2]);
    betaValues  = linspace<vec>(betaInterval[0],betaInterval[1], betaInterval[2]);

}



