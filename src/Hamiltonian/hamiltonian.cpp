#include <src/Hamiltonian/hamiltonian.h>
#include <src/electronInteraction/coulombinteraction.h>
#include <src/electronInteraction/nointeraction.h>

Hamiltonian::Hamiltonian(Kinetic* kinetic,Potential* potential,
                         ElectronInteraction* electronInteraction):
    kinetic(kinetic),
    potential(potential),
    electronInteraction(electronInteraction)
{
}


