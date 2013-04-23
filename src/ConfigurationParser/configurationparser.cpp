#include "configurationparser.h"
#include <src/Minimizer/bfminimizer.h>
#include <src/Minimizer/steepestdescent.h>

ConfigurationParser::ConfigurationParser(const int& nProcess, const int& myRank):
    nProcess(nProcess),
    myRank(myRank)

{
}

/************************************************************
Name:
Description:
*/
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

/************************************************************
Name:
Description:
*/
void ConfigurationParser::singleRun()
{
    vmcapp= new VMCApp(myRank,nProcess);
    vmcapp->loadConfiguration(&cfg);
    vmcapp->alpha=alpha;
    vmcapp->beta=beta;
    vmcapp->runVMCApp(nCycles,idum);
    vmcapp->messagePassing();
}

/************************************************************
Name:
Description:
*/
void ConfigurationParser::chooseAndRunMinimization()
{
    minimizerType = cfg.lookup("setup.MinimizerSettings.minimizerType");

    switch (minimizerType) {
    case BRUTEFORCE:
        minimizer= new BFMinimizer(myRank,nProcess);
        minimizer->loadConfiguration(&cfg);
        minimizer->runMinimizer();
        break;

    case STEEPESTDESCENT:
        minimizer= new SteepestDescent(myRank,nProcess);
        minimizer->loadConfiguration(&cfg);
        minimizer->runMinimizer();
        break;
    }
}

/************************************************************
Name:
Description:
*/
void ConfigurationParser::runBlocking()
{
    if (myRank==0){
        cout << "Starting blocking analysis." << endl;
        Blocking block(nProcess);
        block.loadConfiguration(&cfg);
        block.doBlocking();

        string dataPath = "../vmc/results/blocking";
        // Plotting result
        string pythonPath = "python " + dataPath
                + "/plotBlocking.py "
                + dataPath + "/blocking.mat";
        system(pythonPath.c_str());
    }
}


/************************************************************
Name:               loadConfiguration
Description:        loads different variables
*/
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
    nCycles = cfg.lookup("AppSettings.cycles");
    idum = cfg.lookup("AppSettings.idum");

}


#if ONEBODYDENSITY
    if (myRank==0){
        cout << "Onebody Density" << endl;
    }
    OnebodyDensity onebodyDensity(nProcess, myRank);
    onebodyDensity.loadConfiguration(&cfg);
    onebodyDensity.computeOnebodyDensity();

    if (myRank==0){
        string dataPath = "../vmc/results/onebodyDensity";
        // Plotting result
        string pythonPath = "python " + dataPath
                + "/OBD.py "
                + dataPath + "/OBD.mat";
        system(pythonPath.c_str());
    }
#endif
