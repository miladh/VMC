#include <iostream>
#include <src/VMCApp/vmcapp.h>
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

    VMCApp parser(nProcess, myRank);
    parser.setup();

    MPI_Finalize();

    end = clock();
    timeSpent = (double)(end - begin) / CLOCKS_PER_SEC;
    if(myRank==0){
        cout << "Execution time: " << timeSpent << endl;
    }


    return 0;
}

