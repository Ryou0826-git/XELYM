//gene_dist.cpp

#include "MT.h"
#include <cmath>

double Uniform() {
    return genrand_real3();
}

double rand_exp(double lambda){
    return -log(Uniform())/lambda;
}

double rand_normal(double mu, double sigma){
    double z=sqrt(-2.0*log(Uniform())) * sin(2.0*M_PI*Uniform() );
    return mu + sigma*z;
}

double rand_chi(int k){
    int i;
    double z, w=0;

    for(i = 0; i < k; ++i) {
        z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
        w+=z*z;
    }

    return w;
}

double rand_cauchy(double mu, double gamma){
    return mu + gamma*tan(M_PI*(Uniform()-0.5));
}

double rand_Lnormal(double mu, double sigma){
    double z= mu + sigma*sqrt(-2.0*log(Uniform()))*sin(2.0*M_PI*Uniform());
    return exp(z);
}

double rand_Igauss(double mu, double lambda){
   double x,y,w,z;
   x=sqrt(-2.0*log(Uniform()))*sin(2.0*M_PI*Uniform());  //normal random number
   y=x*x;                                                //chi-squared
   w= mu + 0.5*y*mu*mu/lambda -(0.5*mu/lambda)*sqrt(4.0*mu*lambda*y+mu*mu*y*y);
   z=Uniform();

   if(z< mu/(mu+w)) {
       return w;
   } else {
       return mu*mu/w;
   }
}

double rand_gamma(double theta, double kappa){

    int int_kappa;
    double frac_kappa;

    int_kappa  = (int)kappa;
    frac_kappa = kappa - (double)int_kappa;

    double u, uu;
    double b, p, x_frac, x_int;
    int    i;

    /*integer part*/
    x_int=0;
    for(i=0;i<int_kappa;i++){
        x_int+=-log(Uniform()); // add expnential random number with mean 1
    }

    /*fractional part*/
    if (fabs(frac_kappa) < 0.01) {
        x_frac=0;
    } else {
        b=(exp(1.0)+frac_kappa)/exp(1.0);
        while(1){

            u=Uniform();
            p=b*u;

            uu=Uniform();

            if (p<=1.0) {
                x_frac=pow(p,1.0/frac_kappa);
                if (uu<=exp(-x_frac)) {
                    break;
                }
            }

            else{
                x_frac=-log((b-p)/frac_kappa);
                if (uu<=pow(x_frac,frac_kappa-1.0)) {
                    break;
                }
            }

        }
    }

    return (x_int+x_frac)*theta;
}
