#include <iostream>
#include <libconfig.h++>
#include "src/Minimizer/minimizer.h"
#include "src/slater/slater.h"
#include "src/Blocking/blocking.h"

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"



using namespace std;
using namespace libconfig;

int main()
{
    clock_t begin, end;
    double timeSpent;
    begin = clock();

    Config cfg;
    cfg.readFile("../vmc/src/config.cfg");


    int nProcess, myRank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    Minimizer min(myRank,nProcess);
    min.loadConfiguration(&cfg);
    min.runMinimizaer();

    MPI_Finalize();


    end = clock();
    timeSpent = (double)(end - begin) / CLOCKS_PER_SEC;

    if(myRank==0){
    cout << "Execution time: " <<timeSpent << endl;
    }



#if BLOCKING
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
#endif

    return 0;
}

