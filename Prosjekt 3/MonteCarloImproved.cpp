#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#include "../libraries/lib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "ran0.cpp"
using namespace  std;
// output file as global variable

ofstream ofile;
// Function declaration
double Improved_Monte_Carlo(double *x);\
double gaussian_deviate(long * idum);

int main()
{
    int n = 10000000;
    double x[6], fx ;
    double int_mc = 0.; double variance = 0.;
    double sum_sigma= 0. ; long idum=-1 ;
    double jacobidet=4*pow(acos(-1.),4.)*1./16;
    for ( int i = 1; i <= n; i++){
        for (int j = 0; j< 2; j++) {
            double y=ran0(&idum);
            x[j]=-0.25*log(1.-y);//r coordinates
        }
        for (int j=2;j<4;j++){
            x[j]=2*acos(-1.)*ran0(&idum);//fi coordinates
        }
        for (int j=4;j<6;j++){
            x[j]=acos(-1.)*ran0(&idum);//theta coordinates
        }
        //double y=ran0(&idum);
        fx=Improved_Monte_Carlo(x);
        int_mc += fx;//4*exp(-4.*y);
        sum_sigma += fx*fx;
    }

    cout << setiosflags(ios::showpoint | ios::uppercase);
    int_mc = int_mc/((double)n);
    sum_sigma =  sum_sigma/((double) n);
    variance=sum_sigma-int_mc*int_mc;
    cout << "int_mc  " << jacobidet*int_mc << endl;
    cout << "Sigma  " << jacobidet*sqrt(variance/(double)n) << endl;
    cout << "Variance  " << jacobidet*jacobidet*(variance/(double)n) << endl;
    //return 0;
//    for (int j=0;j<6;j++){

//     cout << x[j] << endl;
//    }
}

double Improved_Monte_Carlo(double *x)
{
//    for (int j=0;j<6;j++){

//     cout << x[j] << endl;
    //}

     double r12=((x[0]*x[0])+(x[1]*x[1])-2*x[0]*x[1]*(cos(x[4])*cos(x[5]) +
            sin(x[4])*sin(x[5])*cos(x[2]-x[3])));
     if (r12 < 1e-8){
         return 0;
     }

     double g = (x[0]*x[0])*(x[1]*x[1])*
                 sin(x[4])*sin(x[5])
                 /sqrt(r12);
     return g;
}


/*double gaussian_deviate(long * idum){
    static int iset=0;
    static double gset;
    double fac, rsq, v1, v2;

    if (idum<0) iset=0;
    if (iset == 0){
        do {
            v1=2.*ran0(idum) -1.0;
            v2=2.*ran0(idum) -1.0;
            rsq = v1*v1+v2*v2;
        }
        while (rsq >= 1.0 || rsq == 0.);
        fac = sqrt(-2.*log(rsq)/rsq);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    }else{
        iset = 0;
        return gset;
    }
}
*/
