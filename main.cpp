
#include "global.h"
#include "parameters.h"
#include "omega.h"
#include "TDMA.h"
#include "solvediffeq.h"
#include "phi.h"
#include "Q_partition.h"
#include "Conc.h"
#include "Incomp.h"
#include "pin.h"
#include "output.h"
#include "fE.h"
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
    int *Ns;
    double *drz;
    double **chiMatrix;
    double fE_hom;
    
    //Allocate memory
    w=create_3d_double_array(ChainType,Nr,Nz,"w");          //Auxiliary potential fields
    eta=create_2d_double_array(Nr,Nz,"eta");                //Incompressibility field
    sigma=create_2d_double_array(Nr, Nz, "sigma");          //Pinning condition field
    phi=create_3d_double_array(ChainType,Nr,Nz,"phi");      //Concentration fields
    chi=create_1d_double_array(ChainType,"chi");            //Interaction parameters
    f=create_1d_double_array(ChainType,"f");                //Chain fractions
    Ns=create_1d_integer_array(ChainType, "Ns");            //Chain lengths
    drz=create_1d_double_array(2,"drz");                    //Grid spacing
    mu=create_1d_double_array(2, "mu");                     //Chemical potentials
    tip=create_1d_integer_array(2,"tip");                   //Pore tip location
    chiMatrix=create_2d_double_array(ChainType,ChainType,"chiMatrix");
    
    //Initial time for random number generator
    long iseed;
    time_t t;
    iseed=time(&t);
    srand48(iseed);
    
    //Set parameters
    parameters(chi,f,ds,Ns,drz,chiMatrix,tip,mu);
    
    //Calculate homogeneous free energy
    fE_hom=homogfE(mu,chiMatrix,f);
    
    //Set up initial omega field
    omega(w);
    
    //SCFT
    FreeEnergy(w,phi,eta,Ns,ds,chi,drz,chiMatrix,sigma,tip,mu,fE_hom);
    
    //Destroy memory allocations------------
    destroy_3d_double_array(w);
    destroy_2d_double_array(eta);
    destroy_3d_double_array(phi);
    destroy_1d_double_array(chi);
    destroy_1d_integer_array(Ns);
    destroy_1d_integer_array(tip);
    destroy_1d_double_array(f);
    destroy_1d_double_array(drz);
    destroy_2d_double_array(chiMatrix);
    //-------------------------------------
    
    return 0;
}
