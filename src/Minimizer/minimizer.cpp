#include <src/Minimizer/minimizer.h>


Minimizer::Minimizer(Config *cfg, const int& myRank, const int& nProcess):
    cfg(cfg),
    myRank(myRank),
    nProcess(nProcess),
    parser(new ConfigurationParser(cfg,myRank,nProcess))

{
}


