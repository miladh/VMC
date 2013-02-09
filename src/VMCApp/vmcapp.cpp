#include "vmcapp.h"
#include "src/VMCSolver/mcbf.h"
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
    charge(2),
    stepLength(1.0),
    h(0.001),
    h2(1000000),
    alpha(1.85),
    beta(0.25),
    idum(-1)


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


}

void VMCApp::run()
{
    vmcsolver = new MCBF();
    vmcsolver->solve(nDimensions,nParticles,hamiltonian,TrialWaveFunction,nCycles,idum);
    double energySquared= vmcsolver->energySquared;
    double energy =vmcsolver->energy;
    cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
}



