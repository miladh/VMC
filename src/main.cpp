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

    MPI_Init(NULL, NULL);

    Minimizer min;
    min.loadConfiguration(&cfg);
    min.runMinimizaer();

    MPI_Finalize();
    return 0;
}

