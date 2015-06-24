
double q_partition(double ***qdagB,double ***qC, double *drz, int *Ns, double *mu){
    
    double Q,Q_AB,Q_C;
    int i,j;
    
    Q=0.0;
    Q_AB=0.0;
    Q_C=0.0;
    
    //corners
    Q_AB+=0.25*qdagB[0][0][Ns[1]]*drz[0]*drz[1]*r_0;
    Q_C+=0.25*qC[0][0][Ns[2]]*drz[0]*drz[1]*r_0;
    
    Q_AB+=0.25*qdagB[(int)Nr-1][0][Ns[1]]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    Q_C+=0.25*qC[(int)Nr-1][0][Ns[2]]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    
    Q_AB+=0.25*qdagB[(int)Nr-1][(int)Nz-1][Ns[1]]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    Q_C+=0.25*qC[(int)Nr-1][(int)Nz-1][Ns[2]]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    
    Q_AB+=0.25*qdagB[0][(int)Nz-1][Ns[1]]*drz[0]*drz[1]*r_0;
    Q_C+=0.25*qC[0][(int)Nz-1][Ns[2]]*drz[0]*drz[1]*r_0;
   
   
    
    //sides
    i=0;
    for (j=1;j<(int)Nz-1;j++){
        Q_AB+=0.5*qdagB[i][j][Ns[1]]*drz[0]*drz[1]*r_0;
        Q_C+=0.5*qC[i][j][Ns[2]]*drz[0]*drz[1]*r_0;
    }
    i=(int)Nr-1;
    for (j=1;j<(int)Nz-1;j++){
        Q_AB+=0.5*qdagB[i][j][Ns[1]]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
        Q_C+=0.5*qC[i][j][Ns[2]]*drz[0]*drz[1]*(drz[0]*((double)Nr-1.0)+r_0);
    }
    j=0;
    for (i=1;i<(int)Nr-1;i++){
        Q_AB+=0.5*qdagB[i][j][Ns[1]]*drz[0]*drz[1]*(drz[0]*((double)i)+r_0);
        Q_C+=0.5*qC[i][j][Ns[2]]*drz[0]*drz[1]*(drz[0]*((double)i)+r_0);
    }
    j=int(Nz-1);
    for (i=1;i<(int)Nr-1;i++){
        Q_AB+=0.5*qdagB[i][j][Ns[1]]*drz[0]*drz[1]*(drz[0]*((double)i)+r_0);
        Q_C+=0.5*qC[i][j][Ns[2]]*drz[0]*drz[1]*(drz[0]*((double)i)+r_0);
    }
    
    //Middle
    for(i=1;i<(int)Nr-1;i++){
        for(j=1;j<(int)Nz-1;j++){
            Q_AB+=qdagB[i][j][Ns[1]]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]);
            Q_C+=qC[i][j][Ns[2]]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]);
        }
    }
    Q_AB=exp(mu[0])*Q_AB;
    Q_C=(exp(mu[1]*kappa)*Q_C)/kappa;
    cout<<"Q_AB: "<<Q_AB<<" Q_C: "<<Q_C<<endl;
    //I'm adding the two single chain partition functions together for the return function
    Q=Q_AB+Q_C;
    
    
    return Q;
}