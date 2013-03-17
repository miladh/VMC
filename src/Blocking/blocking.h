#ifndef BLOCKING_H
#define BLOCKING_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>

using namespace std;
using namespace arma;
using namespace libconfig;

class Blocking
{
public:
    Blocking(const int &nProcess);
    void loadConfiguration(Config *cfg);
    void doBlocking();

protected:
    vec readDataFromFile();
    vec computeBlock(const vec &data, int blockSize);
    string dataPath;
    int deltaBlockSize;
    int nNodes;
    double maxBlockSizeTreshold;
    int stepLength;
    string outFilename;
    string dataName;
};

#endif // BLOCKING_H
