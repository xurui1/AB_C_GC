double homogfE(double *mu, double **chimatrix, double *f){

    int i,j;
    double pA_ave,pB_ave,pC_ave;
    double wA_ave,wB_ave,wC_ave;
    double dwA_ave,dwB_ave,dwC_ave,dpp_ave;
    double eta_ave;
    double f_int, f_omeg,fE_hom;
    double *p_vect;
    double *w_vect;

    p_vect=create_1d_double_array(3,"p_vect");
    w_vect=create_1d_double_array(3,"w_vect");

    f_int=0.0;
    f_omeg=0.0;

    dwA_ave=0.0;
    dwB_ave=0.0;
    dwC_ave=0.0;


    eta_ave=0.0;

    pA_ave=0.002;
    pB_ave=pA_ave;
    pC_ave=1.0-(pA_ave+pB_ave);

    wA_ave=chimatrix[0][1]*pB_ave+chimatrix[0][2]*pC_ave+eta_ave;
    wB_ave=chimatrix[1][0]*pA_ave+chimatrix[1][2]*pC_ave+eta_ave;
    wC_ave=chimatrix[2][0]*pA_ave+chimatrix[2][1]*pB_ave+eta_ave;


    for (i=0;i<10000000;i++){

        eta_ave=eta_ave-0.05*(1.0-(pA_ave+pB_ave+pC_ave));

        pA_ave=exp(mu[0]-wA_ave*f[0]-wB_ave*f[1])*f[0];
        pB_ave=exp(mu[0]-wA_ave*f[0]-wB_ave*f[1])*f[1];
        pC_ave=exp(kappa*(mu[1]-wC_ave));

        dwA_ave=(chimatrix[0][1]*pB_ave+chimatrix[0][2]*pC_ave+eta_ave)-wA_ave;
        dwB_ave=(chimatrix[1][0]*pA_ave+chimatrix[1][2]*pC_ave+eta_ave)-wB_ave;
        dwC_ave=(chimatrix[2][0]*pA_ave+chimatrix[2][1]*pB_ave+eta_ave)-wC_ave;

        dpp_ave=1.0-(pA_ave+pB_ave+pC_ave);

        wA_ave=wA_ave+0.005*dwA_ave;
        wB_ave=wB_ave+0.005*dwB_ave;
        wC_ave=wC_ave+0.005*dwC_ave;

    }
//not sure if I will need these
    //phiAB_hom=pA_ave+pB_ave;
    //phiC_hom=pC_ave;


    p_vect[0]=pA_ave;
    p_vect[1]=pB_ave;
    p_vect[2]=pC_ave;


    w_vect[0]=wA_ave;
    w_vect[1]=wB_ave;
    w_vect[2]=wC_ave;


    for (i=0;i<3;i++){
        for (j=i;j<3;j++){
            f_int+=p_vect[i]*p_vect[j]*chimatrix[i][j];
        }
        f_omeg+=p_vect[i]*w_vect[i];
    }


    fE_hom=f_int-f_omeg-(exp(mu[0]-wA_ave*f[0]-wB_ave*f[1]))-(exp(kappa*(mu[1]-wC_ave))/kappa);

    return fE_hom;
}
