#include "vmcapp.h"
#include "src/Solver/mcbf.h"
#include "src/includes/lib.h"
#include "src/Wavefunction/jastrowwavefunction.h"
#include "src/Wavefunction/basicwavefunction.h"
#include "src/Wavefunction/hydrogenicwavefunction.h"
#include "src/Potential/coulomb_potential.h"
#include "src/Kinetic/numericalkinetic.h"
#include "src/Kinetic/closedformkinetic.h"
//#include <mpi.h>

VMCApp::VMCApp(Config *cfg)
{
    this->cfg=cfg;
}




/************************************************************
Name:               runVMCApp
Description:        starts VMC calculations
*/
void VMCApp::runVMCApp()
{



//    // MPI init
//    int numproc, my_rank;
//    MPI_Init(NULL, NULL);
//    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
//    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

//    idum = idum - my_rank - time(NULL);
//    nCycles /= numproc;
//    double tot_energy,tot_energy_sq;

    TrialWaveFunction = new JastrowWavefunction;
    TrialWaveFunction->alpha=alpha;
    TrialWaveFunction->beta=beta;

    potential = new CoulombPotential(cfg);

    kinetic= new NumericalKinetic(cfg);
    kinetic->wf = TrialWaveFunction;
    kinetic->alpha=alpha;
    kinetic->beta=beta;

    hamiltonian =new Hamiltonian();
    hamiltonian->potential=potential;
    hamiltonian->kinetic=kinetic;

    solver = new MCBF(hamiltonian,TrialWaveFunction);
    solver->loadConfiguration(cfg);
    solver->solve();


//    energy =solver->energy;
//    MPI_Allreduce(&energy, &tot_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    tot_energy /= numproc;

//    energySquared= solver->energySquared;
//    MPI_Allreduce(&energySquared, &tot_energy_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    tot_energy_sq /= numproc;


//    // Printing progress.
//    if (my_rank == 0) {
//        cout << alpha << ", " << beta << " Energy = " << tot_energy
//                << ", Variance = " << tot_energy_sq - tot_energy * tot_energy
//                << ", Sigma = " << sqrt(tot_energy_sq - tot_energy * tot_energy)
//                << ", MC cycles = " << nCycles * numproc
//                << "\n";
//    }
    energy= solver->energy;
    energySquared =solver->energySquared;
}



