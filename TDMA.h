void TDMA(double *x, const size_t N, const double a[], const double b[], const double c[]){
    
    size_t in;
    
    double *cprime;
    cprime=create_1d_double_array(N, "cprime");
    
    cprime[0]= c[0]/b[0];
    x[0]=x[0]/b[0];
    
    for (in=1;in<N;in++){
        double m = 1.0/(b[in]-a[in]*cprime[in]);
        cprime[in] = c[in]*m;
        x[in] = (x[in]-a[in]*x[in-1])*m;
    }
    
    for (in=N-1;in-->0;){
        x[in] = x[in] - cprime[in]*x[in+1];
    }
    
    destroy_1d_double_array(cprime);
}