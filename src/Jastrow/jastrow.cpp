#include "jastrow.h"

Jastrow::Jastrow()
{
}

/************************************************************
Name:               setaValues
Description:        Initiating a matrix with all the spin dependant a-values.
*/
void Jastrow::setaValues(const uint &nParticles){
    a = zeros(nParticles,nParticles);


    for (uint i = 0; i < nParticles; i++) {
        for (uint j = i+1; j <nParticles; j++) {
            if ((i < nParticles / 2 && j >= nParticles / 2) || (i >=nParticles / 2 && j < nParticles / 2)){
                a(i, j) = 0.5;
                a(j,i) = 0.5;
            }
            else{
                a(i, j) = 0.25;
                a(j,i) = 0.25;
        }
        }
    }


//    for (uint i = 0; i < nParticles; i++) {
//        for (uint j = i; j <nParticles; j++) {
//            if (i == j)
//                a(i, j) = 0;
//            else if ((i < nParticles / 2 && j >= nParticles / 2) || (i >=nParticles / 2 && j < nParticles / 2))
//                a(i, j) = 1.0/2.0;
//            else
//                a(i, j) = 1.0 / 4.0;
//        }
//    }
}


/************************************************************
Name:               evaluateJastrow
Description:        computes jastrowfactor
*/
double Jastrow::evaluateJastrow(const mat &r){

    correlation=0;
    for (uint i=0; i<r.n_rows; i++) {
        for (uint j=i+1; j<r.n_rows; j++) {
            rij= norm( r.row(i)-r.row(j) ,2);

            correlation+=a(i,j)*rij/(1+beta*rij);
        }
    }

    JastrowFactor = exp(correlation);
    return JastrowFactor;

}


/************************************************************
Name:               gradientJastrowEvaluate
Description:        Computes the total Jasrow Wavefunction's
                    gradient in r, component i.
*/
rowvec Jastrow::gradientJastrowEvaluate(const mat &r, uint i) {
    double r_ki,b_ij;
    dJastrowFactor = zeros(1,r.n_cols);

    // Before i
    for (uint k = 0; k < i; k++) {
        r_ki = norm(r.row(k) - r.row(i), 2);
        b_ij=(1 + beta * r_ki);
        dJastrowFactor += (a(k, i) / (b_ij*b_ij))*((r.row(i) - r.row(k))/r_ki);
    }

    // After i
    for (uint k = i + 1; k < r.n_rows; k++) {
        r_ki = norm(r.row(k) - r.row(i), 2);
        b_ij=(1 + beta * r_ki);
         dJastrowFactor += (a(k, i) / (b_ij*b_ij))*((r.row(i) - r.row(k))/r_ki);
    }
    return dJastrowFactor;
}


/************************************************************
Name:               laplaceJastrowEvaluate
Description:
*/

double Jastrow::laplaceJastrowEvaluate(const mat &r){

double r_ki,b_ij;
    ddJastrowFactor = 0;


    for(uint i=0; i<r.n_rows; i++){
        // Before i
        for (uint k = 0; k < i; k++) {
            r_ki = norm(r.row(k) - r.row(i), 2);
            b_ij=(1 + beta * r_ki);
            ddJastrowFactor  += 2*a(k, i)/(b_ij*b_ij*b_ij*r_ki);

        }
        // After i
        for (uint k = i + 1; k < r.n_rows; k++) {
            r_ki = norm(r.row(k) - r.row(i), 2);
            b_ij=(1 + beta * r_ki);
            ddJastrowFactor  += 2*a(k, i)/(b_ij*b_ij*b_ij*r_ki);
        }

        ddJastrowFactor += dot(gradientJastrowEvaluate(r,i),gradientJastrowEvaluate(r,i));
    }

    return ddJastrowFactor;
}