#include <iostream>
#include <libconfig.h++>
#include "src/Minimizer/minimizer.h"
#include "src/slater/slater.h"
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

//    Slater slater(2);

//    rowvec r1(1,3), r2(1,3);
//    mat r(2,3);

//    r1 << 1 << 0 <<0;
//    r2 << 1 << 0 <<0;

//    r.row(0)=r1;
//    r.row(1)=r2;
//    cout << slater.evaluateSlater(r) <<endl;


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
    cout << "Execution time: "<<timeSpent << endl;

    return 0;
}

