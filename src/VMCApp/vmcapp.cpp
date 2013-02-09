#include "vmcapp.h"
#include "src/Solver/mcbf.h"
#include "src/includes/lib.h"
#include "src/Wavefunction/jastrowwavefunction.h"
#include "src/Wavefunction/basicwavefunction.h"
#include "src/Potential/coulomb_potential.h"
#include "src/Kinetic/numericalkinetic.h"
#include "src/Kinetic/closedformkinetic.h"


VMCApp::VMCApp():
    nDimensions(3),
    nParticles(2),
    nCycles(1000000),
    idum(-1)

{
}

void VMCApp::runVMCApp()
{
    TrialWaveFunction = new JastrowWaveFunction();
    TrialWaveFunction->alpha=alpha;
    TrialWaveFunction->beta=beta;

    potential = new CoulombPotential(); //This should be dropped if Close-form is used

    kinetic= new ClosedFormKinetic();
    kinetic->wf = TrialWaveFunction;
    kinetic->alpha=alpha;
    kinetic->beta=beta;

    hamiltonian =new Hamiltonian();
    hamiltonian->potential=potential;
    hamiltonian->kinetic=kinetic;

    solver = new MCBF();
    solver->solve(nDimensions,nParticles,hamiltonian,TrialWaveFunction,nCycles,idum);
    energySquared= solver->energySquared;
    energy =solver->energy;
}



