void parameters(double *chi,double *f,double &ds,int *Ns,double *drz,double **chiMatrix,int *tip, double *mu){
    
    int Ds=100;
    
    initial=0;
    
    //Length ratio of c homopolymer to diblock copolymer
    kappa=1.0;
    
    r_0=10.0;   //Distance from i=0 to the centre of the cylindre
    
    //Interaction parameters
    chi[0]=20.0;        //Chi_AB
    chi[1]=20.0;        //Chi_BC
    chi[2]=0.0;         //Chi_AC
    
     //Pore Location
    tip[0]=33;
    tip[1]=24;
    
    //Chemical potential array
    mu[0]=0.0;
    mu[1]=-20.0;
    
    //Chain fraction array
    f[0]=0.5;
    f[1]=1.0-f[0];
    f[2]=kappa*1.0;
    
    //Chain length array
    Ns[0]=(Ds*f[0]);
    Ns[1]=(Ds*f[1]);
    Ns[2]=(Ds*f[2]);
    
    cout<<Ns[0]<<" "<<Ns[1]<<" "<<Ns[2]<<endl;
    
    //Step size in r,z direction
    drz[0]=0.16;
    drz[1]=0.16;
    
    //Step length along polymer
    ds=1.0/Ds;
    
    //Interaction Matrix
    chiMatrix[0][0]=0.0;
    chiMatrix[0][1]=chi[0];
    chiMatrix[0][2]=chi[2];
    chiMatrix[1][0]=chi[0];
    chiMatrix[1][1]=0.0;
    chiMatrix[1][2]=chi[1];
    chiMatrix[2][0]=chi[2];
    chiMatrix[2][1]=chi[1];
    chiMatrix[2][2]=0.0;
    
};
