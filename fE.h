double fE(double ***newW, double ***phi, double **chiMatrix, double *drz){
    
    //Needed to use a trapezoidal integration method - pretty messy
    int i,j,ii,jj;
    double fEW,fEchi;
    double fE_int;
    double volume;
    
    //corners
    i=0;
    j=0;
    for(ii=0;ii<ChainType;ii++){
        for(jj=0;jj<ChainType;jj++){
            fEchi+=0.25*phi[ii][i][j]*chiMatrix[ii][jj]*phi[jj][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]);
        }
        fEW+=0.25*(newW[ii][i][j]*phi[ii][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]));
    }
    
    i=0;
    j=Nz-1;
    for(ii=0;ii<ChainType;ii++){
        for(jj=0;jj<ChainType;jj++){
            fEchi+=0.25*phi[ii][i][j]*chiMatrix[ii][jj]*phi[jj][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]);
        }
        fEW+=0.25*(newW[ii][i][j]*phi[ii][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]));
    }
    i=Nr-1;
    j=Nz-1;
    for(ii=0;ii<ChainType;ii++){
        for(jj=0;jj<ChainType;jj++){
            fEchi+=0.25*phi[ii][i][j]*chiMatrix[ii][jj]*phi[jj][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]);
        }
        fEW+=0.25*(newW[ii][i][j]*phi[ii][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]));
    }
    i=Nr-1;
    j=0;
    for(ii=0;ii<ChainType;ii++){
        for(jj=0;jj<ChainType;jj++){
            fEchi+=0.25*phi[ii][i][j]*chiMatrix[ii][jj]*phi[jj][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]);
        }
        fEW+=0.25*(newW[ii][i][j]*phi[ii][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]));
    }
    
    //sides
    i=0;
    for (j=1;j<Nz-1;j++){
        for(ii=0;ii<ChainType;ii++){
            for(jj=0;jj<ChainType;jj++){
                fEchi+=0.5*phi[ii][i][j]*chiMatrix[ii][jj]*phi[jj][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]);
            }
            fEW+=0.5*(newW[ii][i][j]*phi[ii][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]));
        }
    }
    
    i=Nr-1;
    for (j=1;j<Nz-1;j++){
        for(ii=0;ii<ChainType;ii++){
            for(jj=0;jj<ChainType;jj++){
                fEchi+=0.5*phi[ii][i][j]*chiMatrix[ii][jj]*phi[jj][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]);
            }
            fEW+=0.5*(newW[ii][i][j]*phi[ii][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]));
        }
    }
    
    j=0;
    for (i=1;i<Nr-1;i++){
        for(ii=0;ii<ChainType;ii++){
            for(jj=0;jj<ChainType;jj++){
                fEchi+=0.5*phi[ii][i][j]*chiMatrix[ii][jj]*phi[jj][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]);
            }
            fEW+=0.5*(newW[ii][i][j]*phi[ii][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]));
        }
    }
    j=Nz-1;
    for (i=1;i<Nr-1;i++){
        for(ii=0;ii<ChainType;ii++){
            for(jj=0;jj<ChainType;jj++){
                fEchi+=0.5*phi[ii][i][j]*chiMatrix[ii][jj]*phi[jj][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]);
            }
            fEW+=0.5*(newW[ii][i][j]*phi[ii][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]));
        }
    }
    
    //middle
    for(i=1;i<Nr-1;i++){
        for(j=1;j<Nz-1;j++){
            for(ii=0;ii<ChainType;ii++){
                for(jj=0;jj<ChainType;jj++){
                    fEchi+=phi[ii][i][j]*chiMatrix[ii][jj]*phi[jj][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]);
                }
                fEW+=(newW[ii][i][j]*phi[ii][i][j]*((double)i*drz[0]+r_0)*(drz[0])*(drz[1]));
            }
        }
    }
    //normalize by box size
    volume=Pi*(pow((drz[0]*((double)Nr-1.0)+r_0),2.0)-pow(r_0,2.0))*(drz[1]*((double)Nz-1));
    fEchi/=volume/(Pi);
    fEW/=volume/(2.0*Pi);
    
    fE_int=fEchi-fEW;

    return fE_int;
    
}