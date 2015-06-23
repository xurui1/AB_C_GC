
//Build LHS matrix (non-zero elements) for ADI step #2
void Matrix_r(int ii, double **w, double *drz, double ds, double *rmid, double *rupper, double *rlower){
    int i;
    
    for (i=0;i<Nr;i++){
        rmid[i]=0.0;
        rmid[i]=1.0+(ds/(pow((double)drz[0],2.0)))+(ds/4.0)*w[i][ii];
    }
    
    for (i=1;i<(int)Nr-1;i++){
        rupper[i]=0.0;
        rupper[i]=-(ds/(2.0*pow((double)drz[0],2.0)))-(ds/((((double)drz[0]*((double)i+1.0))+r_0)*4.0*(double)drz[0]));
    }
    for (i=0;i<(int)Nr-2;i++){
        rlower[i]=0.0;
        rlower[i]=-(ds/(2.0*pow((double)drz[0],2.0)))+(ds/((((double)drz[0]*((double)i-1.0))+r_0)*4.0*(double)drz[0]));
    }
    
    rupper[0]=-ds/(pow((double)drz[0],2.0))-(ds/((((double)drz[0]*(1.0))+r_0)*4.0*(double)drz[0]))+(ds/((((double)drz[0]*(-1.0))+r_0)*4.0*(double)drz[0]));
;

    rlower[Nr-2]=-ds/(pow((double)drz[0],2.0))-(ds/((((double)drz[0]*(((double)Nr-1))+r_0)*4.0*(double)drz[0])))+(ds/((((double)drz[0]*((double)Nr-3.0))+r_0)*4.0*(double)drz[0]));

}

//Build LHS matrix (non-zero elements) for ADI step #1
void Matrix_z(int ii, double **w, double *drz, double ds, double *zmid, double *zupper, double *zlower){
    int i;
    
    for (i=0;i<Nz;i++){
        zmid[i]=0.0;
        zmid[i]=1.0+(ds/pow((double)drz[1],2.0))+(ds/4.0)*w[ii][i];
    }
    
    for (i=0;i<(int)Nz-1;i++){
        zupper[i]=0.0;
        zlower[i]=0.0;
        zupper[i]=-(ds/(2.0*pow((double)drz[1],2.0)));
        zlower[i]=-(ds/(2.0*pow((double)drz[1],2.0)));
    }
    zupper[0]=2.0*zupper[0];
    zlower[Nz-2]=2.0*zlower[Nz-2];
    
}

//Apply ADI method to solve the modified diffusion equation in 2D cylindrical coordinates
void solvediffyQ(double ***q, double **w, double **qint, double ds, int Ns, double *drz){

    int            i,j,s;  // some counters
    double gamma, betaL, betaU;
    double *bvecr;
    double *bvecz;
    double *zmid;
    double *zupper;
    double *zlower;
    double *rmid;
    double *rupper;
    double *rlower;
    
    //Allocate memory for `matrix' operations
    bvecr=create_1d_double_array(Nr, "bvecr");
    bvecz=create_1d_double_array(Nz, "bvecz");
    zmid=create_1d_double_array(Nz, "zmid");
    zupper=create_1d_double_array(Nz-1, "zupper");
    zlower=create_1d_double_array(Nz-1, "zlower");
    rmid=create_1d_double_array(Nr, "rmid");
    rupper=create_1d_double_array(Nr-1, "rupper");
    rlower=create_1d_double_array(Nr-1, "rlower");
    
    
    for (s=0;s<(int)Ns+1;s++){
        
        //Empty RHS vector
        for (j=0;j<Nz;j++){
            bvecz[j]=0.0;
        }
        
        //Step 1-scan over z for all r
        for (i=0;i<Nr;i++){
            Matrix_z(i,w,drz,ds,zmid,zupper,zlower);
            
            //Create b-vector (bvecz) for step 1
            for (j=0;j<Nz;j++){
                gamma=1.0-(ds/(pow((double)drz[0],2.0)))-((ds/4.0)*w[i][j]);
                betaL=(ds/(2.0*pow((double)drz[0],2.0)))-(ds/(((((double)i-1.0)*(double)drz[0])+(r_0))*4.0*(double)drz[0]));
                betaU=(ds/(2.0*pow((double)drz[0],2.0)))+(ds/(((((double)i+1.0)*(double)drz[0])+(r_0))*4.0*(double)drz[0]));
                //cout<<i<<" "<<j<<" "<<s<<" g: "<<gamma<<" bL: "<<betaL<<" bU: "<<betaU<<endl;
                if (i==0){
                    bvecz[j]=gamma*qint[i][j]+betaL*qint[i+1][j]+betaU*qint[i+1][j];
                }
                else if (i==(int)Nr-1){
                    bvecz[j]=gamma*qint[i][j]+betaL*qint[i-1][j]+betaU*qint[i-1][j];
                }
                else{
                    bvecz[j]=gamma*qint[i][j]+betaL*qint[i-1][j]+betaU*qint[i+1][j];
                }
            }
            //Use TDMA to solve matrix algebra problem
            TDMA(bvecz,Nz,zlower,zmid,zupper);
            
            for (j=0;j<Nz;j++){
                q[i][j][s]=bvecz[j];
                //cout<<i<<" "<<j<<" "<<s<<" bvecz: "<<bvecz[j]<<endl;
                bvecz[j]=0.0;
            }
            
        }
        
        
        //Intermediate Step
        for (i=0;i<Nr;i++){
            for (j=0;j<Nz;j++){
                qint[i][j]=q[i][j][s];
                //cout<<"i: "<<i<<" j: "<<j<<" s: "<<s<<" qint: "<<q[i][j][s]<<endl;
            }
        }
        
        //Empty RHS vector #2
        for (i=0;i<Nr;i++){
            bvecr[i]=0.0;
        }
        
        //Step 2-scan over i for all z
        for (j=0;j<Nz;j++){
            Matrix_r(j,w,drz,ds,rmid,rupper,rlower);
            
               for (i=0;i<Nr;i++){
                   gamma=1.0-(ds/(pow((double)drz[1],2.0)))-((ds/4.0)*w[i][j]);
                   betaL=ds/(2.0*pow((double)drz[1],2.0));
                   betaU=ds/(2.0*pow((double)drz[1],2.0));
                   if(j==0){
                       bvecr[i]=gamma*qint[i][j]+betaL*qint[i][j+1]+betaU*qint[i][j+1];
                   }
                   else if(j==(int)Nz-1){
                       bvecr[i]=gamma*qint[i][j]+betaL*qint[i][j-1]+betaU*qint[i][j-1];
                   }
                   else{
                       bvecr[i]=gamma*qint[i][j]+betaL*qint[i][j-1]+betaU*qint[i][j+1];
                   }
               }
            //Use TDMA to solve matrix algebra problem
            TDMA(bvecr,Nr,rlower,rmid,rupper);
            
            //Now we have our solution for all i,j for s, from s-1. Full step completed.
            for (i=0;i<Nr;i++){
                q[i][j][s]=bvecr[i];
                //cout<<"i: "<<i<<" j: "<<j<<" s: "<<s<<" q: "<<q[i][j][s]<<endl;
                qint[i][j]=bvecr[i];
                bvecr[i]=0.0;
            }
        }

    }
    
    
    //Deallocate memory
    destroy_1d_double_array(bvecr);
    destroy_1d_double_array(bvecz);
    destroy_1d_double_array(zmid);
    destroy_1d_double_array(zupper);
    destroy_1d_double_array(zlower);
    destroy_1d_double_array(rmid);
    destroy_1d_double_array(rlower);
    destroy_1d_double_array(rupper);


}

