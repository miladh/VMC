#include <iostream>
#include <libconfig.h++>
#include <src/Minimizer/bfminimizer.h>
#include <src/Minimizer/steepestdescent.h>
#include <src/Blocking/blocking.h>
#include <src/OnebodyDensity/onebodydensity.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"



using namespace std;
using namespace libconfig;

int main()
{

    int nProcess, myRank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    clock_t begin, end;


    double timeSpent;
    begin = clock();

    Config cfg;

    for(int i=0; i<nProcess;i++){
        if(myRank==i){
            cfg.readFile("../vmc/src/config.cfg");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


//    BFMinimizer min(myRank,nProcess);
//    min.loadConfiguration(&cfg);
//    min.runMinimizer();
    SteepestDescent min(myRank,nProcess);
    min.loadConfiguration(&cfg);
    min.runMinimizer();



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



#if BLOCKING
    if (myRank==0){
        cout << "Starting blocking analysis." << endl;
    }
    Blocking block(nProcess);
    block.loadConfiguration(&cfg);
    block.doBlocking();|

    if (myRank==0){
        string dataPath = "../vmc/results/blocking";
        // Plotting result
        string pythonPath = "python " + dataPath
                + "/plotBlocking.py "
                + dataPath + "/blocking.mat";
        system(pythonPath.c_str());
    }
#endif

    MPI_Finalize();

    end = clock();
    timeSpent = (double)(end - begin) / CLOCKS_PER_SEC;
    if(myRank==0){
        cout << "Execution time: " << timeSpent << endl;
    }

    return 0;
}

