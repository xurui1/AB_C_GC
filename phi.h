void phi_total(double **phiA, double **phiB, double **phiC, double *drz){
    
    //Here I am calculating the total concentration of each species using a trapezoidal (?) rule. I removed the trapezoidal component because I don't think it should matter too much if we are using the Neumann boundary condition. Alternatively, I think I would have to modify the incompressibility condition.
    
    double phiA_tot,phiB_tot,phiC_tot;
    double volume;
    phiA_tot=0.0;
    phiB_tot=0.0;
    phiC_tot=0.0;
    double phi_tot=0.0;
    int i,j;
    
    
    //corners
    phiA_tot+=0.25*phiA[0][0]*drz[0]*drz[1]*r_0;
    phiB_tot+=0.25*phiB[0][0]*drz[0]*drz[1]*r_0;
    phiC_tot+=0.25*phiC[0][0]*drz[0]*drz[1]*r_0;
    
    phiA_tot+=0.25*phiA[(int)Nr-1][0]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    phiB_tot+=0.25*phiB[(int)Nr-1][0]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    phiC_tot+=0.25*phiC[(int)Nr-1][0]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    
    phiA_tot+=0.25*phiA[(int)Nr-1][(int)Nz-1]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    phiB_tot+=0.25*phiB[(int)Nr-1][(int)Nz-1]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    phiC_tot+=0.25*phiC[(int)Nr-1][(int)Nz-1]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    
    phiA_tot+=0.25*phiA[0][(int)Nz-1]*drz[0]*drz[1]*r_0;
    phiB_tot+=0.25*phiB[0][(int)Nz-1]*drz[0]*drz[1]*r_0;
    phiC_tot+=0.25*phiC[0][(int)Nz-1]*drz[0]*drz[1]*r_0;
    
    //sides
    i=0;
    for (j=1;j<(int)Nz-1;j++){
        phiA_tot+=0.5*phiA[i][j]*drz[0]*drz[1]*r_0;
        phiB_tot+=0.5*phiB[i][j]*drz[0]*drz[1]*r_0;
        phiC_tot+=0.5*phiC[i][j]*drz[0]*drz[1]*r_0;
    }
    i=(int)Nr-1;
    for (j=1;j<(int)Nz-1;j++){
        phiA_tot+=0.5*phiA[i][j]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
        phiB_tot+=0.5*phiB[i][j]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
        phiC_tot+=0.5*phiC[i][j]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    }
    j=0;
    for (i=1;i<(int)Nr-1;i++){
        phiA_tot+=0.5*phiA[i][j]*drz[0]*drz[1]*(drz[0]*((double)i)+r_0);
        phiB_tot+=0.5*phiB[i][j]*drz[0]*drz[1]*(drz[0]*((double)i)+r_0);
        phiC_tot+=0.5*phiC[i][j]*drz[0]*drz[1]*(drz[0]*((double)i)+r_0);
    }
    j=int(Nz-1);
    for (i=1;i<(int)Nr-1;i++){
        phiA_tot+=0.5*phiA[i][j]*drz[0]*drz[1]*(drz[0]*((double)i)+r_0);
        phiB_tot+=0.5*phiB[i][j]*drz[0]*drz[1]*(drz[0]*((double)i)+r_0);
        phiC_tot+=0.5*phiC[i][j]*drz[0]*drz[1]*(drz[0]*((double)i)+r_0);
    }
    
    //centre
    
    //This is the way that the incompressibility condition is calculated, although I'm possibly making a mistake here.
    
    for (i=1;i<(int)Nr-1;i++){
        for (j=1;j<(int)Nz-1;j++){
            phiA_tot+=drz[1]*(drz[0]*(double)i+r_0)*drz[0]*phiA[i][j];
            phiB_tot+=drz[1]*(drz[0]*(double)i+r_0)*drz[0]*phiB[i][j];
            phiC_tot+=drz[1]*(drz[0]*(double)i+r_0)*drz[0]*phiC[i][j];
        }
    }
    
    //Volume calculation
    volume=Pi*(pow((drz[0]*((double)Nr-1.0)+r_0),2.0)-pow(r_0,2.0))*(drz[1]*((double)Nz-1));
    //normalize
    phiA_tot/=volume/(2.0*Pi);
    phiB_tot/=volume/(2.0*Pi);
    phiC_tot/=volume/(2.0*Pi);
    phi_tot=phiA_tot+phiB_tot+phiC_tot;
    
    cout<<"phiA: "<<phiA_tot<<" phiB: "<<phiB_tot<<" phiC: "<<phiC_tot<<" total: "<<phi_tot<<endl;
    
    
}


void phi_calc(double ***phi,double ***qA,double ***qdagA,double ***qB,double ***qdagB,double ***qC,int *Ns,double*mu,double ds){
    int i,j,s;
    
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            
            //Empty array elements
            phi[0][i][j]=0.0;
            phi[1][i][j]=0.0;
            phi[2][i][j]=0.0;
            
            //phiA integration
            for(s=0;s<(int)Ns[0]+1;s++){
                if(s==0 || s==(int)Ns[0]){
                    phi[0][i][j]+=0.5*qA[i][j][s]*qdagA[i][j][Ns[0]-s]*ds;
                }
                else{
                    phi[0][i][j]+=qA[i][j][s]*qdagA[i][j][Ns[0]-s]*ds;
                }
            }
            
            //phiB integration
            for(s=0;s<(int)Ns[1]+1;s++){
                if(s==0 || s==(int)Ns[1]){
                    phi[1][i][j]+=0.5*qB[i][j][s]*qdagB[i][j][Ns[1]-s]*ds;
                }
                else{
                    phi[1][i][j]+=qB[i][j][s]*qdagB[i][j][Ns[1]-s]*ds;
                }
            }
            
            //phiC integration
            for(s=0;s<(int)Ns[2]+1;s++){
                if(s==0 || s==(int)Ns[2]){
                    phi[2][i][j]+=0.5*qC[i][j][s]*qC[i][j][Ns[2]-s]*ds;
                }
                else{
                    phi[2][i][j]+=qC[i][j][s]*qC[i][j][Ns[2]-s]*ds;
                }
            }
            
            //Grand canonical relation
            phi[0][i][j]=exp(mu[0])*phi[0][i][j];
            phi[1][i][j]=exp(mu[0])*phi[1][i][j];
            phi[2][i][j]=exp((mu[1])*kappa)*phi[2][i][j]*(1.0/kappa);
        }
    }
    
}
