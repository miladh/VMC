#include <src/Minimizer/steepestdescent.h>

SteepestDescent::SteepestDescent(Config *cfg, const int &myRank, const int &nProcess):
    Minimizer(cfg,myRank,nProcess)
{
    loadAndSetConfiguration();
}

//*****************************************************************************
void SteepestDescent::runMinimizer()
{

    myfile.open ("../vmc/DATA/results");
    myfile << "Alpha   " << "Beta    " << "Energy       " << "Variance    "
           << "Sigma       "<< "Acceptance    "<< endl;


    if(diatomicSystem){
        Setting& RVector = cfg->lookup("setup.MinimizerSettings.BFMinSettings.R");
        vec RInterval = zeros(3);

        for(int i =0; i < RVector.getLength(); i++){
            RInterval[i] = RVector[i];
        }

        vec RValues = linspace<vec>(RInterval[0],RInterval[1],RInterval[2]);
        parameters.resize(nVariationalParameters+1);


        for(uint i = 0; i < RValues.n_elem; i++){
            parameters[nVariationalParameters] = RValues[i];
            minimize();
        }

    }else{
        minimize();
    }

    myfile.close();
}


//*****************************************************************************
void SteepestDescent::minimize()
{
    vec variationalDerivate    = zeros<vec>(nVariationalParameters);
    vec variationalDerivateOld = zeros<vec>(nVariationalParameters);
    vec parametersOld          = zeros<vec>(nVariationalParameters);

    for(int j = 0; j < maxIteration; j++){

        parser->setVariationalParameters(parameters);
        parser->runSolver();

        variationalDerivate = parser->getVariationalDerivate();


        for(int i = 0; i < nVariationalParameters; i++){

            parametersOld[i] = parameters[i];

            if(variationalDerivateOld(i)/variationalDerivate(i) >= 0){
                step *=1.25;
                parameters[i] -= step*signFunc(variationalDerivate(i));
            }else{
                step*=0.5;
                parameters[i] -= step*signFunc(variationalDerivate(i));
            }

        }

        if( fabs(parametersOld[0] - parameters[0]) < epsilon &&
                fabs(parametersOld[1] - parameters[1]) < epsilon){

            if(myRank == 0){
                getResultsAndWrite();
                cout << "Optimal parameters: "<< parametersOld[0]
                     << "  " << parametersOld[1] << endl;
            }

            break;
        }

        variationalDerivateOld = variationalDerivate;
    }

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
void SteepestDescent::loadAndSetConfiguration()
{
    diatomicSystem = cfg->lookup("setup.system");
    step  = cfg->lookup("setup.MinimizerSettings.SDMinSettings.step");
    maxIteration =cfg->lookup("setup.MinimizerSettings.SDMinSettings.maxIteration");
    epsilon = cfg->lookup("setup.MinimizerSettings.SDMinSettings.epsilon");
    Setting& variationalParamters = cfg->lookup("setup.MinimizerSettings.SDMinSettings.varParameters");
    nVariationalParameters = variationalParamters.getLength();

    for(int i = 0; i < nVariationalParameters ; i++){
        parameters.push_back(variationalParamters[i]);
    }


}
