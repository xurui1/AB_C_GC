double Conc(double ***phi,double ***w,int *Ns,double ds,double *drz, double *mu){
    
    int         i,j,s;
    double      Q,Q_AB,Q_C;
    double      volume;
    double      ***qA,***qB,***qC,***qdagA,***qdagB;
    double      **qintA,**qintB,**qintC;
    
    
    
    //Forwards propagators
    qA=create_3d_double_array(Nr,Nz,Ns[0]+1,"qA");
    qB=create_3d_double_array(Nr,Nz,Ns[1]+1,"qB");
    qC=create_3d_double_array(Nr,Nz,Ns[2]+1,"qC");
    
    //Complementary propagators
    qdagA=create_3d_double_array(Nr,Nz,Ns[0]+1,"qdagA");
    qdagB=create_3d_double_array(Nr,Nz,Ns[1]+1,"qdagB");
    
    
    //Intermediate steps for ADI
    qintA=create_2d_double_array(Nr,Nz,"qintA");
    qintB=create_2d_double_array(Nr,Nz,"qintB");
    qintC=create_2d_double_array(Nr,Nz,"qintC");
    
    
    // Here is the for loop for doing the initial conditions setting it to 1.0
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            qintA[i][j]=1.0;
            qintB[i][j]=1.0;
            qintC[i][j]=1.0;
            //cout<<"A: "<<qintA[i][j]<<endl;
        }
    }
    
    
    
    // Here we solve the diffusion equation for the forwards propagators
    solvediffyQ(qA,w[0],qintA,ds,Ns[0],drz);
    solvediffyQ(qB,w[1],qintB,ds,Ns[1],drz);
    solvediffyQ(qC,w[2],qintC,ds,Ns[2],drz);
    
    // The result from the above calculation becomes qdags initial cond
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
                qintA[i][j]=qB[i][j][Ns[1]];
                qdagA[i][j][0]=qB[i][j][Ns[1]];
                qintB[i][j]=qA[i][j][Ns[0]];
                qdagB[i][j][0]=qA[i][j][Ns[0]];
                //std::cout<<qintA[i][j][l]<< "----"<<qintB[i][j][l] <<std::endl;
        }
    }
    
    // Here we will solve the diffusion equation for the complementory qs
    solvediffyQ(qdagA,w[0],qintA,ds,Ns[0],drz);
    solvediffyQ(qdagB,w[1],qintB,ds,Ns[1],drz);
    
    // Here we are doing the sum to get the single chain partition function
    Q=0.0;
    Q_AB=0.0;
    Q_C=0.0;
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            Q_AB+=qdagB[i][j][Ns[1]]*((double)i*(double)drz[0]+r_0)*((double)drz[0])*((double)drz[1]);
            Q_C+=qC[i][j][Ns[2]]*((double)i*(double)drz[0]+r_0)*((double)drz[0])*((double)drz[1]);
        }
    }
    Q_AB=exp(mu[0])*Q_AB;
    Q_C=(exp(mu[1]*kappa)*Q_C)/kappa;
    cout<<"Q_AB: "<<Q_AB<<" Q_C: "<<Q_C<<endl;
    //I'm adding the two single chain partition functions together for the return function
    Q=Q_AB+Q_C;

    // Normalizing with respect to the volume of the box
    volume=(pow((drz[0]*((double)Nr-1.0)+r_0),2.0)-pow(r_0,2.0))*(drz[1]*(double)Nz);
    Q/=volume;
    
    cout<<"Q: "<< Q<<endl;
    
    // Here we do the concentration calculation by integration over box and chain
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            
                //Empty array elements
                phi[0][i][j]=0.0;
                phi[1][i][j]=0.0;
                phi[2][i][j]=0.0;
            
                //phiA integration
                for(s=0;s<(int)(Ns[0]+1);s++){
                    if(s==0 || s==(int)Ns[0]){
                        phi[0][i][j]+=0.5*qA[i][j][s]*qdagA[i][j][Ns[0]-s]*ds;
                    }else{
                        phi[0][i][j]+=qA[i][j][s]*qdagA[i][j][Ns[0]-s]*ds;
                    }
                }
            
                //phiB integration
                for(s=0;s<(int)(Ns[1]+1);s++){
                    if(s==0 || s==(int)Ns[1]){
                        phi[1][i][j]+=0.5*qB[i][j][s]*qdagB[i][j][Ns[1]-s]*ds;
                    }else{
                        phi[1][i][j]+=qB[i][j][s]*qdagB[i][j][Ns[1]-s]*ds;
                    }
                }
            
                //phiC integration
                for(s=0;s<(int)(Ns[2]+1);s++){
                    if(s==0 || s==(int)Ns[2]){
                        phi[2][i][j]+=0.5*qC[i][j][s]*qC[i][j][Ns[2]-s]*ds;
                    }else{
                        phi[2][i][j]+=qC[i][j][s]*qC[i][j][Ns[2]-s]*ds;
                    }
                }
            
                //Grand canonical relation
                phi[0][i][j]=exp(mu[0])*phi[0][i][j];
                phi[1][i][j]=exp(mu[0])*phi[1][i][j];
                phi[2][i][j]=exp(mu[1]*kappa)*phi[2][i][j]*(1.0/kappa);
        }
    }
    
    //calculation of average concentrations over entire computation box
    phi_calc(phi[0],phi[1],phi[2],drz);
    
    
    
    //clearing the memory
    destroy_3d_double_array(qA);
    destroy_3d_double_array(qB);
    destroy_3d_double_array(qC);
    destroy_3d_double_array(qdagA);
    destroy_3d_double_array(qdagB);
    destroy_2d_double_array(qintA);
    destroy_2d_double_array(qintB);
    destroy_2d_double_array(qintC);
    
    return Q;
    
};