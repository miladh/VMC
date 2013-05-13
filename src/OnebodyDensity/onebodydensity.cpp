#include "onebodydensity.h"
#include "src/includes/lib.h"
#include "src/includes/Defines.h"
#include "src/Jastrow/padejastrow.h"
#include "src/Jastrow/nojastrow.h"

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
#pragma GCC diagnostic warning "-Wunused-parameter"


OnebodyDensity::OnebodyDensity(const int &nProcess, const int &myRank):
    nNodes(nProcess),
    myRank(myRank)
{
}

/************************************************************
Name:
Description:
*/
void OnebodyDensity::computeOnebodyDensity()
{
    idum = idum - myRank - time(NULL);
    r=zeros(nParticles,nDimensions);
    rho=zeros(nSteps);

//    wf=setWavefunction();
    nCycles/=nNodes;

    for (int i = 0; i < nSteps; i++) {
        rho(i) = McIntegrator();
        r(0,0)+=dr;
    }

    if (myRank==0){
        normalize();
        writeToFile();
    }
}


/************************************************************
Name:
Description:
*/
void OnebodyDensity::normalize(){

    double norm = sum(rho);
    rho/=norm;

}


/************************************************************
Name:
Description:
*/
double OnebodyDensity::McIntegrator(){

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

        wfValue=wf->evaluateWavefunction(r)*r(0,0);
        rho+=wfValue*wfValue;
    }


    double tmp = rho;
    MPI_Allreduce(&tmp, &rho, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    rho/= nNodes/nCycles;
    return rho;
}

/************************************************************
Name:
Description:
*/
void OnebodyDensity::writeToFile() {

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
void OnebodyDensity::setWavefunction(){

//    Wavefunction* wf = new Wavefunction();

//    switch (WavefunctionType) {
//    case JASTROW:
//        wf = new HLikeWavefunction(nParticles);
//        wf->jas=new PadeJastrow(nParticles,beta);
//        wf->jas->setaValues(nParticles);
//        wf->slater->orbitals->k=alpha;
//        break;

//    case BASIC:
//        wf = new HLikeWavefunction(nParticles);
//        wf->jas=new NoJastrow;
//        wf->slater->orbitals->k=alpha;
//        break;

////    case  HYDROGENIC:
////        wf = new HLikeWavefunction(nParticles);
////        wf->jas=new NoJastrow(nParticles);
////        wf->slater->orbitals->k=charge;
////        break;
//    }

//    return wf;
}

/************************************************************
Name:
Description:
*/
void OnebodyDensity::loadConfiguration(Config *cfg){
    nParticles  = cfg->lookup("setup.nParticles");
    nDimensions = cfg->lookup("setup.nDimensions");
    nCycles=cfg->lookup("OnebodyDensitySettrings.cycles");
    WavefunctionType= cfg->lookup("AppSettings.wavefunction");
    alpha=cfg->lookup("MinimizerSettings.minalpha");
    beta=cfg->lookup("MinimizerSettings.minbeta");
    idum=cfg->lookup("AppSettings.idum");
    charge=cfg->lookup("PotentialSettings.charge");
    nSteps= cfg->lookup("OnebodyDensitySettrings.nSteps");
    dr=cfg->lookup("OnebodyDensitySettrings.dr");
    a=cfg->lookup("OnebodyDensitySettrings.a");
    b=cfg->lookup("OnebodyDensitySettrings.b");

}














