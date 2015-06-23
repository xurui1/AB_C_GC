void phi_calc(double **phiA, double **phiB, double **phiC, double *drz){
    
    //Here I am calculating the total concentration of each species using a trapezoidal (?) rule. I removed the trapezoidal component because I don't think it should matter too much if we are using the Neumann boundary condition. Alternatively, I think I would have to modify the incompressibility condition.
    
    double phiA_tot,phiB_tot,phiC_tot;
    double volume;
    phiA_tot=0.0;
    phiB_tot=0.0;
    phiC_tot=0.0;
    double phi_tot=0.0;
    int i,j;
    
    /*
    //corners
    phiA_tot+=0.25*phiA[0][0];
    phiB_tot+=0.25*phiB[0][0];
    phiC_tot+=0.25*phiC[0][0];
    
    phiA_tot+=0.25*phiA[int(Nr-1)][0];
    phiB_tot+=0.25*phiB[int(Nr-1)][0];
    phiC_tot+=0.25*phiC[int(Nr-1)][0];
    
    phiA_tot+=0.25*phiA[int(Nr-1)][int(Nz-1)];
    phiB_tot+=0.25*phiB[int(Nr-1)][int(Nz-1)];
    phiC_tot+=0.25*phiC[int(Nr-1)][int(Nz-1)];
    
    phiA_tot+=0.25*phiA[0][int(Nz-1)];
    phiB_tot+=0.25*phiB[0][int(Nz-1)];
    phiC_tot+=0.25*phiC[0][int(Nz-1)];
    
    //sides
    i=0;
    for (j=1;j<int(Nz-1);j++){
        phiA_tot+=0.5*phiA[i][j];
        phiB_tot+=0.5*phiB[i][j];
        phiC_tot+=0.5*phiC[i][j];
    }
    i=int(Nr-1);
    for (j=1;j<int(Nz-1);j++){
        phiA_tot+=0.5*phiA[i][j];
        phiB_tot+=0.5*phiB[i][j];
        phiC_tot+=0.5*phiC[i][j];
    }
    j=0;
    for (i=1;i<int(Nr-1);i++){
        phiA_tot+=0.5*phiA[i][j];
        phiB_tot+=0.5*phiB[i][j];
        phiC_tot+=0.5*phiC[i][j];
    }
    j=int(Nz-1);
    for (i=1;i<int(Nr-1);i++){
        phiA_tot+=0.5*phiA[i][j];
        phiB_tot+=0.5*phiB[i][j];
        phiC_tot+=0.5*phiC[i][j];
    }
    */
    //centre
    
    //This is the way that the incompressibility condition is calculated, although I'm possibly making a mistake here.
    
    for (i=0;i<Nr;i++){
        for (j=0;j<Nz;j++){
            phiA_tot+=drz[1]*(drz[0]*(double)i+r_0)*drz[0]*phiA[i][j];
            phiB_tot+=drz[1]*(drz[0]*(double)i+r_0)*drz[0]*phiB[i][j];
            phiC_tot+=drz[1]*(drz[0]*(double)i+r_0)*drz[0]*phiC[i][j];
        }
    }
    
    //Volume calculation
    volume=(pow((drz[0]*((double)Nr-1.0)+r_0),2.0)-pow(r_0,2.0))*(drz[1]*(double)Nz);
    //normalize
    phiA_tot/=volume;
    phiB_tot/=volume;
    phiC_tot/=volume;
    phi_tot=phiA_tot+phiB_tot+phiC_tot;
    
    cout<<"phiA: "<<phiA_tot<<" phiB: "<<phiB_tot<<" phiC: "<<phiC_tot<<" total: "<<phi_tot<<endl;
    
    
}