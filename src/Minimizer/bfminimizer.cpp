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

    for(uint i = 0; i < alphaValues.n_elem; i++){
        parameters[0] = alphaValues[i];

        for(uint j = 0; j < betaValues.n_elem; j++){
            parameters[1] = betaValues[j];

            for(uint k = 0; k < RValues.n_elem; k++){
                parameters[2] = RValues[k];

                parser->setVariationalParameters(parameters);
                parser->runSolver();
                getResultsAndWrite();
            }
        }
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
            myfile << parameters[0] <<"     "<<   parameters[0] <<"     "<<energy
               <<"     "<< variance <<"     "<<sigma
               <<"     "<< acceptance<< endl;
        }
}

//*****************************************************************************
void BFMinimizer::loadAndSetConfiguration()
{
    Setting& intervalPar1 = cfg->lookup("setup.MinimizerSettings.BFMinSettings.alpha");
    Setting& intervalPar2 = cfg->lookup("setup.MinimizerSettings.BFMinSettings.beta");
    Setting& intervalPar3 = cfg->lookup("setup.MinimizerSettings.BFMinSettings.R");

    for(uint i =0; true; i++){
        try{
            double value = 0;
            value = intervalPar1[i];
            alpha.push_back(value);

            value = intervalPar2[i];
            beta.push_back(value);

            value = intervalPar3[i];
            R.push_back(value);

        }catch(exception) {
            break;
        }
    }
    parameters.resize(3);
    alphaValues = linspace<vec>(alpha[0],alpha[1],alpha[2]);
    betaValues  = linspace<vec>(beta[0],beta[1], beta[2]);
    RValues     = linspace<vec>(R[0],R[1],R[2]);

}



