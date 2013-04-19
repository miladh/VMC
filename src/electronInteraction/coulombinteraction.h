#ifndef COULOMBINTERACTION_H
#define COULOMBINTERACTION_H

#include <src/electronInteraction/electroninteraction.h>

class CoulombInteraction : public ElectronInteraction
{
public:
    CoulombInteraction();

    double evaluate(const mat& r);

private:
    double rij,interactionEnergy;
};

#endif // COULOMBINTERACTION_H
