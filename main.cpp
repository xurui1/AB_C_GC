
#include "global.h"
#include "parameters.h"
#include "omega.h"
#include "TDMA.h"
#include "solvediffeq.h"
#include "Conc.h"
#include "Incomp.h"
#include "pin.h"
#include "output.h"
#include "FreeEnergy.h"
#include "homogfE.h"



int main( ){
    
    int *tip;
    double ***w;
    double **eta;
    double **sigma;
    double ***phi;
    double *chi;
    double *f;
    double *mu;
    double ds;
    double *Ns;
    double *drz;
    double **chiMatrix;
    double fE_hom;
    
    w=create_3d_double_array(ChainType,Nr,Nz,"w");
    eta=create_2d_double_array(Nr,Nz,"eta");
    sigma=create_2d_double_array(Nr, Nz, "sigma");
    phi=create_3d_double_array(ChainType,Nr,Nz,"phi");
    chi=create_1d_double_array(ChainType,"chi");
    f=create_1d_double_array(ChainType,"f");
    Ns=create_1d_double_array(ChainType,"Ns");
    drz=create_1d_double_array(2,"drz");
    mu=create_1d_double_array(2, "mu");
    tip=create_1d_integer_array(2,"tip");
    chiMatrix=create_2d_double_array(ChainType,ChainType,"chiMatrix");
    
    long iseed;
    time_t t;
    iseed=time(&t);
    srand48(iseed);
    
    parameters(chi,f,ds,Ns,drz,chiMatrix,tip,mu);
    
    fE_hom=homogfE(mu,chiMatrix,f);
    
    omega(w);
    
    FreeEnergy(w,phi,eta,Ns,ds,chi,drz,chiMatrix,sigma,tip,mu,fE_hom);
    
    
    
    //Destroy memory allocations------------
    
    destroy_3d_double_array(w);
    destroy_2d_double_array(eta);
    destroy_3d_double_array(phi);
    destroy_1d_double_array(chi);
    destroy_1d_double_array(Ns);
    destroy_1d_double_array(f);
    destroy_1d_double_array(drz);
    destroy_2d_double_array(chiMatrix);
    //-------------------------------------
    
    
    return 0;
}
