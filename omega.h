void omega(double ***w){
    
    int i,j,x;
    
    //This is for a sinusoidal initial omega field. Adjust 'initial' in the parameters header.
    //Currently will give a lamellar structure in z-direction
    if (initial==0){
    for(i=0;i<Nr;i++){
        for(j=0;j<Nz;j++){
            w[0][i][j]=-5.0*cos(4.0*Pi*(double)j/(double)Nz)+(double)((rand() %100)/100.0);
            w[1][i][j]=5.0*cos(4.0*Pi*(double)j/(double)Nz)+(double)((rand() %100)/100.0);
            w[2][i][j]=5.0+(double)((rand() %100)/100.0);
        }
    }
    }
    
    //This is for a bilayer conformation
    else if (initial==1){
        ifstream Init1;
        Init1.open("bilayer_M50_N50.dat");
        for (i=0;i<Nr;i++){
            for (j=0;j<Nz;j++){
                Init1 >> x >> w[0][i][j]>>w[1][i][j]>>w[2][i][j];
            }
        }
        Init1.close();
    }
    
    //This is for a disk conformation
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
    
    //This is for a conformation based on the omega field currently stored in results folder
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
    
    //This is for a pore conformation - change tip location if using this one after the disk.
    else if (initial==4){
        ifstream Init4;
        Init4.open("pore_M50_N50.dat");
        for (i=0;i<Nr;i++){
            for (j=0;j<Nz;j++){
                Init4 >> x >> w[0][i][j]>>w[1][i][j]>>w[2][i][j];
            }
        }
        Init4.close();
    }
    
    for (i=0;i<Nr;i++){
            for (j=0;j<Nz;j++){
                std::cout << w[0][i][j]<<" "<<w[1][i][j]<<" "<<w[2][i][j]<<endl;
            }
    }
    
};