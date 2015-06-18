void FreeEnergy(double ***w, double ***phi, double **eta, double *Ns, double ds, double *chi, double *drz, double **chiMatrix, double **sigma, int *tip, double *mu,double fE_hom){
    
    
    double  currentfE, oldfE, deltafE;
    int     maxIter=10000;
    double precision=1e-6;
    int     i,j,iter,chain,ii,jj;
    double  Q;
    double  fEW, fEchi, fES;
    double  epsilon, gamma;
    double  **delphi;
    double  ***delW;
    double  ***newW;
    double  deltaW;
    
    delW=create_3d_double_array(ChainType,Nr,Nz,"delW");
    delphi=create_2d_double_array(Nr,Nz,"delphi");
    newW=create_3d_double_array(ChainType,Nr,Nz,"newW");
    
    currentfE=0.0;
    deltafE=0.0;
    
    epsilon=0.01;
    gamma=0.01;
    
    iter=0;
    std::ofstream outputFile("./results/fE.dat");
    do{
        iter++;
        
        fEW=0.0;
        fEchi=0.0;
        fES=0.0;
        
        Q=Conc(phi,w,Ns,ds,drz,mu);
        
        Incomp(eta,phi,delphi);
        Pin(sigma,phi,tip);
        output(drz,phi,w);
        
        
        fEW=0.0;
        fEchi=0.0;
        deltaW=0.0;

        
        for(i=0;i<Nr;i++){
            for(j=0;j<Nz;j++){
                    for(ii=0;ii<ChainType;ii++){
                        newW[ii][i][j]=0.0;
                }
            }
        }
        
        for(i=0;i<Nr;i++){
            for(j=0;j<Nz;j++){
                
                    for(ii=0;ii<ChainType;ii++){
                        for(jj=0;jj<ChainType;jj++){
                            
                            newW[ii][i][j]+=((chiMatrix[ii][jj]*phi[jj][i][j])+eta[i][j]+sigma[i][j]);
                            fEchi+=phi[ii][i][j]*chiMatrix[ii][jj]*phi[jj][i][j]*drz[0]*drz[1];
                        }
                        fEW+=(newW[ii][i][j]*phi[ii][i][j]*drz[0]*drz[1]);
                        delW[ii][i][j]=newW[ii][i][j]-w[ii][i][j];
                        deltaW+=fabs(delW[ii][i][j]);
                }
            }
        }
        deltaW/=(Nr*Nz);
        fEchi/=(2.0*((Nr*drz[0])*(Nz*drz[1])));
        fEW/=(((Nr*drz[0])*(Nz*drz[1])));
        
        fES=Q;
        
        oldfE=currentfE;
        currentfE=-fES-fEW+fEchi;
        
        deltafE=fabs(currentfE-oldfE);
        
        std::cout<<iter<<" fE:"<<currentfE<< " dfE:"<<currentfE-fE_hom<<" " << deltaW<<std::endl;
        outputFile << iter << " " << currentfE<< " " << currentfE-fE_hom<<" "<< deltaW<<std::endl;
        
        for(i=0;i<Nr;i++){
            for(j=0;j<Nz;j++){
                    for(chain=0;chain<ChainType;chain++){
                        w[chain][i][j]+=(gamma*delW[chain][i][j]-epsilon*delphi[i][j]);
                }
            }
        }
        
        
    }while((iter<maxIter)&&(deltafE>precision)&&(deltaW>precision));
    
    outputFile.close();
    
    destroy_2d_double_array(delphi);
    destroy_3d_double_array(delW);
    destroy_3d_double_array(newW);
    
};
