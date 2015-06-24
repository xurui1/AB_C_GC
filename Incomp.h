void Incomp(double **eta, double ***phi, double **delphi){
    
    int     i,j;
    int     chain;
    double  ptot;
    
    ptot=0.0;
    
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
        
                ptot=0.0;
                delphi[i][j]=0.0;
                
                for(chain=0;chain<ChainType;chain++){
                    ptot+=phi[chain][i][j];
                }
                            
                delphi[i][j]=1.0-ptot;
                eta[i][j]-=delphi[i][j];
                
            }
        }
    
    
};
