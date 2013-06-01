#include "vmcapp.h"
#include <src/Minimizer/bfminimizer.h>
#include <src/Minimizer/steepestdescent.h>

VMCApp::VMCApp(const int& nProcess, const int& myRank):
    nProcess(nProcess),
    myRank(myRank)
{
        loadConfiguration();
}

//*****************************************************************************
void VMCApp::options()
{

    int singleRunIsEnable    = cfg.lookup("setup.singleRun");
    int minimizationIsEnable = cfg.lookup("setup.minimization");
    int blockingIsEnable     = cfg.lookup("setup.blocking");

    if(singleRunIsEnable){
        singleRun();
    }
    else if(minimizationIsEnable){
        chooseAndRunMinimization();
    }
    else if(blockingIsEnable){
        runBlocking();
    }
}

//*****************************************************************************
void VMCApp::singleRun()
{

    if (myRank==0){
        cout << "--------Starting single VMC run-----------" << endl;
    }
    ConfigurationParser* parser = new ConfigurationParser(&cfg,myRank,nProcess);
    parser->runSolver();
}

//*****************************************************************************
void VMCApp::chooseAndRunMinimization()
{

    if (myRank==0){
        cout << "----------Starting Minimizer----------" << endl;
    }

    int minimizerType = cfg.lookup("setup.MinimizerSettings.minimizerType");
    Minimizer* minimizer;

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
void VMCApp::runBlocking()
{
    if (myRank==0){
        cout << "----------Starting blocking analysis----------" << endl;
        Blocking block(nProcess);
        block.loadConfiguration(&cfg);
        block.doBlocking();

//        string dataPath = "../vmc/DATA/blocking";
//        // Plotting result
//        string pythonPath = "python " + dataPath
//                + "/plotBlocking.py "
//                + dataPath + "/blocking.mat";
//        system(pythonPath.c_str());
    }
}


//*****************************************************************************
void VMCApp::loadConfiguration()
{
    for(int i=0; i<nProcess;i++){
        if(myRank==i){
            cfg.readFile("../vmc/src/config.cfg");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
