#include <iostream>
#include <libconfig.h++>
#include "src/Minimizer/minimizer.h"

using namespace std;
using namespace libconfig;
int main()
{

    Config cfg;
    cfg.readFile("../vmc/src/config.cfg");

    Minimizer min;
    min.loadConfiguration(&cfg);
    min.runMinimizaer();

    return 0;
}

