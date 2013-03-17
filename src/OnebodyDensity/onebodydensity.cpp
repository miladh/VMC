#include "onebodydensity.h"
#include "src/includes/lib.h"
#include "src/includes/Defines.h"
#include "src/Wavefunction/hlikewavefunction.h"
#include "src/Jastrow/padejastrow.h"
#include "src/Jastrow/nojastrow.h"

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
#pragma GCC diagnostic warning "-Wunused-parameter"


OnebodyDensity::OnebodyDensity(const int &nProcess, const int &myRank):
    nNodes(nProcess),
    myRank(myRank),
    nSteps(100),
    dr(0.05),
    a(-2),
    b(2)
{
}

/************************************************************
Name:
Description:
*/
void OnebodyDensity::computeOnebodyDensity()
{
    mat r = zeros(nParticles,nDimensions);
    vec rho = zeros(nSteps);
    idum = idum - myRank - time(NULL);

    wf=setWavefunction();
    nCycles/=nNodes;

    for (int i = 0; i < nSteps; i++) {
        rho(i) = McIntegrator(r);
        r(0,0)+=dr;
    }

    normalize(rho);
    writeToFile(r,rho);
}


/************************************************************
Name:
Description:
*/
vec OnebodyDensity::normalize(vec rho){

    double norm = sum(rho);
    rho/=norm;

    return rho;
}


/************************************************************
Name:
Description:
*/
double OnebodyDensity::McIntegrator(mat r){

    double rho=0.0;
    //Lopp over MC cycles
    for (int c=0; c<nCycles; c++) {
        //move the other paricles around
        //while the first electron is stationary
        for (uint i=1; i<nParticles; i++) {
            for(uint j=0; j< nDimensions; j++){
                r(i,j)= a+((b-a)*ran2(&idum));
            }
        }

        wfValue=wf->wavefunction(r)*r(0,0)*r(0,0);
        rho+=wfValue*wfValue;
    }

    double tmp = rho;
    MPI_Allreduce(&tmp, &rho, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    rho/= nNodes;
    return rho;
}

/************************************************************
Name:
Description:
*/
void OnebodyDensity::writeToFile(mat r,vec rho) {

    myfile.open("../vmc/results/onebodyDensity/OBD.mat");
    r(0,0)=0.0;

    //Writing result in file
    for (int i =0; i< nSteps; i++ ){
        myfile << r(0,0) << "      " << rho(i)<< endl;
        r(0,0)+= dr;
    }

}


/************************************************************
Name:
Description:
*/
Wavefunction* OnebodyDensity::setWavefunction(){

    Wavefunction* wf;

    switch (WavefunctionType) {
    case Jastrow:
        wf = new HLikeWavefunction(nParticles);
        wf->jas=new PadeJastrow(nParticles);
        wf->jas->alpha=alpha;
        wf->jas->beta=beta;
        wf->jas->setaValues(nParticles);
        wf->slater->orbitals->k=alpha;
        break;

    case Basic:
        wf = new HLikeWavefunction(nParticles);
        wf->jas=new NoJastrow(nParticles);
        wf->slater->orbitals->k=alpha;
        break;

    case  Hydrogenic:
        wf = new HLikeWavefunction(nParticles);
        wf->jas=new NoJastrow(nParticles);
        wf->slater->orbitals->k=charge;
        break;
    }

    return wf;
}

/************************************************************
Name:
Description:
*/
void OnebodyDensity::loadConfiguration(Config *cfg){
    nDimensions=cfg->lookup("SolverSettings.dim");
    nParticles=cfg->lookup("SolverSettings.N");
    nCycles=cfg->lookup("AppSettings.cycles");
    WavefunctionType= cfg->lookup("AppSettings.wavefunction");
    alpha=cfg->lookup("MinimizerSettings.minalpha");
    beta=cfg->lookup("MinimizerSettings.minbeta");
    idum=cfg->lookup("AppSettings.idum");
    charge=cfg->lookup("PotentialSettings.charge");

}














