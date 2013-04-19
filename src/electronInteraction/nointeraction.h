#ifndef NOINTERACTION_H
#define NOINTERACTION_H

#include <src/electronInteraction/electroninteraction.h>

class NoInteraction : public ElectronInteraction
{
public:
    NoInteraction();

    inline double evaluate(const mat& ){return 0;}
};

#endif // NOINTERACTION_H
