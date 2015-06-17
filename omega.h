void omega(double ***w){
    
    int i,j,chain,x;
    
    if (initial==0){
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            for(chain=0;chain<ChainType;chain++){
                w[chain][i][j]=-5.0*cos(2.*Pi*j/Nr);
            }
        }
    }
    }
    
    else if (initial==1){
        ifstream Init1;
        Init1.open("bilayer_M50_N50.dat");
    for (i=0;i<Nr;i++){
        for (j=0;j<Nz;j++){
            Init1 >> x >> w[0][i][j]>>w[1][i][j]>>w[2][i][j];
        }
    }
    Init1.close();
    
    for (i=0;i<Nr;i++){
        for (j=0;j<Nz;j++){
            std::cout << w[0][i][j]<<" "<<w[1][i][j]<<" "<<w[2][i][j]<<endl;
        }
    }
    }
    
    else if (initial==2){
        ifstream Init2;
        Init2.open("disk_M50_N50.dat");
        for (i=0;i<Nr;i++){
            for (j=0;j<Nz;j++){
                Init2 >> x >> w[0][i][j]>>w[1][i][j]>>w[2][i][j];
            }
        }
        
        
        Init2.close();
    }
    
    else if (initial==3){
        ifstream Init3;
        Init3.open("./results/omega.dat");
        for (i=0;i<Nr;i++){
            for (j=0;j<Nz;j++){
                Init3 >> x >> w[0][i][j]>>w[1][i][j]>>w[2][i][j];
            }
        }
        
        
        Init3.close();
        
    }
    
    
    
};