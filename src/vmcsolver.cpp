#include "vmcsolver.h"
#include "includes/lib.h"

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

VMCSolver::VMCSolver() :
    nDimensions(3),
    charge(2),
    stepLength(1.0),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(-1),
    alpha(1.85),
    beta(0.25),
    nCycles(1000000)
{
}

void VMCSolver::runMonteCarloIntegration()
{
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

    double energySum = 0;
    double energySquaredSum = 0;

    double deltaE;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {

        // Store the current value of the wave function
        waveFunctionOld = TrialWaveFunction.waveFunction(nDimensions,nParticles,rOld,alpha,beta);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }

            // Recalculate the value of the wave function
            waveFunctionNew = TrialWaveFunction.waveFunction(nDimensions,nParticles,rNew,alpha,beta);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    waveFunctionOld = waveFunctionNew;
                }
            } else {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
            }
            // update energies
            deltaE = localEnergy(rNew);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
}

double VMCSolver::localEnergy(const mat &r)
{
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = TrialWaveFunction.waveFunction(nDimensions,nParticles,r,alpha,beta);

    // Kinetic energy

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = TrialWaveFunction.waveFunction(nDimensions,nParticles,rMinus,alpha,beta);
            waveFunctionPlus = TrialWaveFunction.waveFunction(nDimensions,nParticles,rPlus,alpha,beta);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            potentialEnergy += 1 / sqrt(r12);
        }
    }

    return kineticEnergy + potentialEnergy;
}

//double VMCSolver::waveFunction(const mat &r)
//{
//    double rParticle;
//    double correlation, argument;
//    double r12;

//    correlation=argument=0.0;


//    // The correlation factor
//    for (int i=0; i<nParticles-1; i++) {
//        for (int j=i+1; j<nParticles; j++) {
//            r12= 0.0;

//            for (int k=0; k <nDimensions ; k++ ){

//                r12 += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
//            }

//            correlation+=sqrt(r12)/(2+2*beta*sqrt(r12));
//        }
//    }

//    for (int i=0; i<nParticles; i++) {
//        rParticle=0;

//        for (int j=0; j<nDimensions; j++) {
//            rParticle += r(i,j)*r(i,j);
//        }
//        argument += sqrt(rParticle);
//    }

//    //Trial wave function
//    double Trial_func = exp(-alpha* argument)*exp(correlation);

//    return Trial_func;

////    double argument = 0;
////    for(int i = 0; i < nParticles; i++) {
////        double rSingleParticle = 0;
////        for(int j = 0; j < nDimensions; j++) {
////            rSingleParticle += r(i,j) * r(i,j);
////        }
////        argument += sqrt(rSingleParticle);
////    }
////    return exp(-argument * alpha);
//}
