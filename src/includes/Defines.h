#ifndef DEFINES_H
#define DEFINES_H

// 0: Blocking off
// 1: Blocking on
#define BLOCKING 0


// 0: OBD off
// 1: OBD on
#define ONEBODYDENSITY 0

enum SYSTEM{
   ATOMS, MOLECULES
};

enum TrialWavefunction {
   JASTROW, BASIC, HYDROGENIC
};

enum SOLVER {
   BF, IS
};

enum INTERACTION {
   NOINTERACTION, COULOMBINTERACTION
};


enum MINIMIZER {
   BRUTEFORCE, STEEPESTDESCENT
};
#endif // DEFINES_H
