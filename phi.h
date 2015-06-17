void phi_calc(double **phiA, double **phiB, double **phiC, double *drz){
    
    //Here I am calculating the total concentration of each species using a trapezoidal (?) rule
    
    double phiA_tot,phiB_tot,phiC_tot;
    phiA_tot=0.0;
    phiB_tot=0.0;
    phiC_tot=0.0;
    double phi_tot;
    int i,j;
    
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
    
    //centre
    
    for (i=1;i<int(Nr-1);i++){
        for (j=1;j<int(Nz-1);j++){
            phiA_tot+=phiA[i][j];
            phiB_tot+=phiB[i][j];
            phiC_tot+=phiC[i][j];
        }
    }
    
    phiA_tot/=Nr*Nz;
    phiB_tot/=Nr*Nz;
    phiC_tot/=Nr*Nz;
    phi_tot=phiA_tot+phiB_tot+phiC_tot;
    
    cout<<"phiA: "<<phiA_tot<<" phiB: "<<phiB_tot<<" phiC: "<<phiC_tot<<" total: "<<phi_tot<<endl;
    
    
}