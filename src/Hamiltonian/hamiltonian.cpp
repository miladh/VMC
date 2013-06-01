#include <src/Hamiltonian/hamiltonian.h>
#include <src/electronInteraction/coulombinteraction.h>
#include <src/electronInteraction/nointeraction.h>

Hamiltonian::Hamiltonian(Config *cfg, Kinetic* kinetic, Potential* potential,
                         ElectronInteraction* electronInteraction):
    cfg(cfg),
    kinetic(kinetic),
    potential(potential),
    electronInteraction(electronInteraction)

{
}


