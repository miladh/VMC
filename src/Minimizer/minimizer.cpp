#include <src/Minimizer/minimizer.h>


Minimizer::Minimizer(Config *cfg, const int& myRank, const int& nProcess):
    cfg(cfg),
    myRank(myRank),
    nProcess(nProcess),
    vmcapp(new VMCApp(cfg,myRank,nProcess))

{
}


