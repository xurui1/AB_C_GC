double Conc(double ***phi,double ***w,double *Ns,double ds,double *drz, double *mu){
    
    int         i,j,s;
    double      Q,Q_AB,Q_C;
    double      ***qA,***qB,***qC,***qdagA,***qdagB;
    double      **qint,**qintA,**qintB,**qintC;
    
    
    
    
    qA=create_3d_double_array(Nr,Nz,((int)Ns[0]+1),"qA");
    qB=create_3d_double_array(Nr,Nz,((int)Ns[1]+1),"qB");
    qC=create_3d_double_array(Nr,Nz,((int)Ns[2]+1),"qC");
    qdagA=create_3d_double_array(Nr,Nz,((int)Ns[0]+1),"qdagA");
    qdagB=create_3d_double_array(Nr,Nz,((int)Ns[1]+1),"qdagB");
    qint=create_2d_double_array(Nr,Nz,"qint");
    qintA=create_2d_double_array(Nr,Nz,"qintA");
    qintB=create_2d_double_array(Nr,Nz,"qintB");
    qintC=create_2d_double_array(Nr,Nz,"qintC");
    
    
    // Here is the for loop for doing the initial conditions setting it to 1.0
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            qintA[i][j]=1.0;
            qintB[i][j]=1.0;
            qintC[i][j]=1.0;
            //qC[i][j][0]=1.0;
        }
    }
    
    // Here we will solve the diffusion question
    solvediffyQ(qA,w[0],qintA,ds,(int)Ns[0],drz);
    solvediffyQ(qB,w[1],qintB,ds,(int)Ns[1],drz);
    solvediffyQ(qC,w[2],qintC,ds,(int)Ns[2],drz);
    
    // The result from the above calculation becomes qdags initial cond
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
                qintA[i][j]=qB[i][j][(int)Ns[0]];
                qintB[i][j]=qA[i][j][(int)Ns[1]];
                //std::cout<<qintA[i][j][l]<< "----"<<qintB[i][j][l] <<std::endl;
        }
    }
    
    // Here we will solve the diffusion equation for the complementory qs
    solvediffyQ(qdagA,w[0],qintA,ds,(int)Ns[0],drz);
    solvediffyQ(qdagB,w[1],qintB,ds,(int)Ns[1],drz);
    //solvediffyQ(qC,w[2],qintC,ds,(int)Ns[2],drz);
    
    // Here we are doing the sum to get the single chain partition function
    Q=0.0;
    Q_AB=0.0;
    Q_C=0.0;
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            Q_AB+=qdagB[i][j][(int)Ns[1]]*drz[0]*drz[1];
            Q_C+=qC[i][j][(int)Ns[2]]*drz[0]*drz[1];
        }
    }
    Q=exp(mu[0])*Q_AB+(exp(mu[1]*kappa)*Q_C)/kappa;
    // Normalizing with respect to the volume of the box
    Q/=((drz[0]*Nr)*(drz[1]*Nz));
    
    
    // Here we do the concentration calculation by integration over box and chain
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            
                phi[0][i][j]=0.0;
                phi[1][i][j]=0.0;
                phi[2][i][j]=0.0;
                
                for(s=0;s<(Ns[0]+1);s++){
                    if(s==0 || s==(int)Ns[0]){
                        phi[0][i][j]+=0.5*qA[i][j][s]*qdagA[i][j][(int)Ns[0]-s]*ds;
                    }else{
                        phi[0][i][j]+=qA[i][j][s]*qdagA[i][j][(int)Ns[0]-s]*ds;
                    }
                }
                for(s=0;s<(Ns[1]+1);s++){
                    if(s==0 || s==(int)Ns[1]){
                        phi[1][i][j]+=0.5*qB[i][j][s]*qdagB[i][j][(int)Ns[1]-s]*ds;
                    }else{
                        phi[1][i][j]+=qB[i][j][s]*qdagB[i][j][(int)Ns[1]-s]*ds;
                    }
                }
            
                for(s=0;s<(Ns[2]+1);s++){
                    if(s==0 || s==(int)Ns[2]){
                        phi[2][i][j]+=0.5*qC[i][j][s]*qC[i][j][(int)Ns[2]-s]*ds;
                    }else{
                        phi[2][i][j]+=qC[i][j][s]*qC[i][j][(int)Ns[2]-s]*ds;
                    }
                }
                phi[0][i][j]=exp(mu[0])*phi[0][i][j];
                phi[1][i][j]=exp(mu[0])*phi[1][i][j];
                phi[2][i][j]=exp(mu[1]*kappa)*phi[2][i][j]*(1.0/kappa);
        }
    }
    
    phi_calc(phi[0],phi[1],phi[2],drz);
    
    
    
    //clearing the memory
    destroy_3d_double_array(qA);
    destroy_3d_double_array(qB);
    destroy_3d_double_array(qC);
    destroy_3d_double_array(qdagA);
    destroy_3d_double_array(qdagB);
    destroy_2d_double_array(qint);
    destroy_2d_double_array(qintA);
    destroy_2d_double_array(qintB);
    destroy_2d_double_array(qintC);
    
    
    return Q;
    
    
};