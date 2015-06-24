void Pin(double **sigma, double ***phi, int *tip){
    
    int i,j,Ntip,Mtip;
    
    for (i=0;i<Nr;i++){
        for (j=0;j<Nz;j++){
            sigma[i][j]=0.0;
        }
    }
    
    if (initial==2 || initial==4){
        Ntip=tip[0];
        Mtip=tip[1];
        sigma[Ntip][Mtip]=sigma[Ntip][Mtip]-10.0*(phi[0][Ntip][Mtip]+phi[2][Ntip][Mtip]-phi[1][Ntip][Mtip]);
    }
}