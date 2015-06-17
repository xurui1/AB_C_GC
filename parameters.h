void parameters(double *chi,double *f,double &ds,double *Ns,double *drz,double **chiMatrix,int *tip, double *mu){
    
    int Ds=100;
    initial=1;
    
    bond_energy=0.2;
    N_bond=100;
    kappa=1.0;
    z1=1.0;
    z2=1.0;
    r_0=10.0;
    
    chi[0]=30.0;
    chi[1]=30.0;
    chi[2]=0.0;
    
    tip[0]=20;
    tip[1]=24;
    
    mu[0]=-1.0;
    mu[1]=0.0;
    
    f[0]=0.5;
    f[1]=1.0-f[0];
    f[2]=kappa*1.0;
    
    Ns[0]=(int)(Ds*f[0]);
    Ns[1]=(int)(Ds*f[1]);
    Ns[2]=(int)(Ds*f[2]);
    
    drz[0]=0.12;
    drz[1]=0.12;
    
    ds=1.0/Ds;
    
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
