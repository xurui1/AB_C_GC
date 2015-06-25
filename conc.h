double Conc(double ***phi,double ***w,int *Ns,double ds,double *drz, double *mu){
    
    int         i,j,s;
    double      Q;
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
    
    
    // Here is the for loop for setting the propagator initial conditions to 1.0
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            qintA[i][j]=1.0;
            qA[i][j][0]=1.0;
           
            qintB[i][j]=1.0;
            qB[i][j][0]=1.0;
            
            qintC[i][j]=1.0;
            qC[i][j][0]=1.0;
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
    
    // Here we are get the single chain partition functions Q_AB+Q_C
    Q=q_partition(qdagB,qC,drz,Ns,mu);
    

    // Normalizing with respect to box volume
    volume=Pi*(pow((drz[0]*((double)Nr-1.0)+r_0),2.0)-pow(r_0,2.0))*(drz[1]*((double)Nz-1.0));
    Q/=volume/(2.0*Pi);
    
    cout<<"Q: "<< Q<<endl;
    
    // Here we do the concentration calculation by integration over box and chain
    phi_calc(phi,qA,qdagA,qB,qdagB,qC,Ns,mu,ds);

    //calculation of average concentrations over entire computation box
    phi_total(phi[0],phi[1],phi[2],drz);
    
    
    
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