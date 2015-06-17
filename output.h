void output(double *drz, double ***phi,double ***w){
    
    int i,j;
    
    ofstream outputFile1("./results/phi.dat");
    
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            outputFile1 <<i*drz[0]<<" "<<j*drz[1]<<" "<<phi[0][i][j]<<" "<<phi[1][i][j]<<" "<<phi[2][i][j]<<std::endl;
        }
    }
    
    outputFile1.close();
    
    ofstream outputFile2("./results/omega.dat");
    
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            outputFile2 <<i*drz[0]<<" "<<j*drz[1]<<" "<<w[0][i][j]<<" "<<w[1][i][j]<<" "<<w[2][i][j]<<std::endl;
        }
    }
    outputFile2.close();

    
}