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
double brute_force_MC(double *x);
double ran0(long *idum);
int main()
{
    int n = 10000000;
    double x[6], y, fx ;
    double int_mc = 0.; double variance = 0.;
    double sum_sigma= 0. ; long idum=-1 ;
    double length=1.5;
    double jacobidet=pow((2*length),6.);
    //evaluate the integral with importance sampling
    for ( int i = 1; i <= n; i++){

        for (int j = 0; j< 6; j++) {

            //x[j]=-length+2*length*ran0(&idum);
            x[j]=-length+2*length*ran0(&idum);
        }
        fx=brute_force_MC(x);
        int_mc += fx;
        sum_sigma += fx*fx;
    }
    int_mc = int_mc/((double)n);
    sum_sigma =  sum_sigma/((double) n );
    variance=sum_sigma-int_mc*int_mc;
    cout << "int_mc  " << jacobidet*int_mc << endl;
    cout << "Sigma  " << jacobidet*sqrt(variance/(double)n) << endl;
    cout << "Variance  " << variance << endl;

    return 0;
}

double brute_force_MC(double *x)

{
    double alpha = 2.;
    // evaluate the different terms of the exponential
    double exp1 = -2*alpha*sqrt((x[0]*x[0])+(x[1]*x[1])+(x[2]*x[2]));
    double exp2 = -2*alpha*sqrt((x[3]*x[3])+(x[4]*x[4])+(x[5]*x[5]));
    double deno = sqrt(pow((x[0]-x[3]),2)+pow((x[1]-x[4]),2)+pow((x[2]-x[5]),2));
    double value = exp(exp1+exp2)/deno;
    return value;

}


//#define IA 16807
//#define IM 2147483647
//#define AM (1.0/IM)
//#define IQ 127773
//#define IR 2836
//#define MASK 123459876

//double ran0(long *idum)
//{
//   long     k;
//   double   ans;

//   *idum ^= MASK;
//   k = (*idum)/IQ;
//   *idum = IA*(*idum - k*IQ) - IR*k;
//   if(*idum < 0) *idum += IM;
//   ans=AM*(*idum);
//   *idum ^= MASK;
//   return ans;
//}
//#undef IA
//#undef IM
//#undef AM
//#undef IQ
//#undef IR
//#undef MASK
