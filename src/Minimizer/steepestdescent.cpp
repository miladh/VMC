#include <src/Minimizer/steepestdescent.h>

SteepestDescent::SteepestDescent(Config *cfg, const int &myRank, const int &nProcess):
    Minimizer(cfg,myRank,nProcess)
{
    loadAndSetConfiguration();
}

//*****************************************************************************
void SteepestDescent::runMinimizer()
{
    vec old = zeros<vec>(parameters.size());

    myfile.open ("../vmc/DATA/results");
    myfile << "Alpha   " << "Beta    " << "Energy       " << "Variance    "
           << "Sigma       "<< "Acceptance    "<< endl;


    for(int i=0; i < 100; i++){

        parser->setVariationalParameters(parameters);
        parser->runSolver();

        variationalDerivate = parser->getVariationalDerivate();
        getResultsAndWrite();


        for(uint i = 0; i < parameters.size(); i++){

            if(old(i)/variationalDerivate(i) >= 0){
                step *=1.25;
                parameters[i] -= step*signFunc(variationalDerivate(i));
            }else{
                step*=0.5;
                parameters[i] -= step*signFunc(variationalDerivate(i));
            }

        }

        old = variationalDerivate;
    }
        myfile.close();
}


//*****************************************************************************
int SteepestDescent::signFunc(double varDer)
{
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
    myfile <<  parameters[0] <<"     "<<   parameters[0] <<"     "<<energy
           <<"     "<< variance <<"     "<<sigma
           <<"     "<< acceptance<< endl;
    }
}

//*****************************************************************************
void SteepestDescent::loadAndSetConfiguration()
{
    step  = cfg->lookup("setup.MinimizerSettings.SDMinSettings.step");
    Setting& variationalParamters = cfg->lookup("setup.MinimizerSettings.SDMinSettings.varParameters");

    for(uint i =0; true; i++){
        try{
            double value = 0;
            value = variationalParamters[i];
            parameters.push_back(value);

        }catch(exception) {
            break;
        }
    }

    variationalDerivate = zeros<vec>(parameters.size());

}
