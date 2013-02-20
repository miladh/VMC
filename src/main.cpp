#include <iostream>
#include <libconfig.h++>
#include "src/Minimizer/minimizer.h"
#include <mpi.h>

using namespace std;
using namespace libconfig;
int main()
{

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
    return 0;
}

