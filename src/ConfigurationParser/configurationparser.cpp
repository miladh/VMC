#include "configurationparser.h"
#include <src/Minimizer/bfminimizer.h>
#include <src/Minimizer/steepestdescent.h>

ConfigurationParser::ConfigurationParser(const int& nProcess, const int& myRank):
    nProcess(nProcess),
    myRank(myRank)
{
}

//*****************************************************************************
void ConfigurationParser::setup()
{
    loadConfiguration();
    if(singleRunIsEnable){
        singleRun();
    }
    if(minimizationIsEnable){
        chooseAndRunMinimization();
    }
    if(blockingIsEnable){
        runBlocking();
    }
}

//*****************************************************************************
void ConfigurationParser::singleRun()
{
    vmcapp = new VMCApp(&cfg,myRank,nProcess);
    vmcapp->alpha = alpha;
    vmcapp->beta  = beta;
    vmcapp->runVMCApp();
}

//*****************************************************************************
void ConfigurationParser::chooseAndRunMinimization()
{
    minimizerType = cfg.lookup("setup.MinimizerSettings.minimizerType");

    switch (minimizerType) {
    case BRUTEFORCE:
        minimizer= new BFMinimizer(&cfg, myRank,nProcess);
        minimizer->runMinimizer();
        break;

    case STEEPESTDESCENT:
        minimizer= new SteepestDescent(&cfg, myRank,nProcess);
        minimizer->runMinimizer();
        break;
    }
}

//*****************************************************************************
void ConfigurationParser::runBlocking()
{
    if (myRank==0){
        cout << "Starting blocking analysis." << endl;
        Blocking block(nProcess);
        block.loadConfiguration(&cfg);
        block.doBlocking();

        string dataPath = "../vmc/DATA/blocking";
        // Plotting result
        string pythonPath = "python " + dataPath
                + "/plotBlocking.py "
                + dataPath + "/blocking.mat";
        system(pythonPath.c_str());
    }
}


//*****************************************************************************
void ConfigurationParser::loadConfiguration()
{
    for(int i=0; i<nProcess;i++){
        if(myRank==i){
            cfg.readFile("../vmc/src/config.cfg");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    singleRunIsEnable = cfg.lookup("setup.singleRun");
    minimizationIsEnable =cfg.lookup("setup.minimization");
    blockingIsEnable = cfg.lookup("setup.blocking");

    alpha =cfg.lookup("setup.singleRunSettings.alpha");
    beta = cfg.lookup("setup.singleRunSettings.beta");

}


#if ONEBODYDENSITY
    if (myRank==0){
        cout << "Onebody Density" << endl;
    }
    OnebodyDensity onebodyDensity(nProcess, myRank);
    onebodyDensity.loadConfiguration(&cfg);
    onebodyDensity.computeOnebodyDensity();

    if (myRank==0){
        string dataPath = "../vmc/DATA/onebodyDensity";
        // Plotting result
        string pythonPath = "python " + dataPath
                + "/OBD.py "
                + dataPath + "/OBD.mat";
        system(pythonPath.c_str());
    }
#endif
