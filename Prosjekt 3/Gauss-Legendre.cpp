#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#include "lib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "gauleg.cpp"
#include "time.h"
using namespace  std;
// output file as global variable
ofstream ofile;

double int_function(double x1, double y1, double z1, double x2, double y2, double z2);
void output(int N, double a, double b, double int_gauss);


int main()
{
    int N = 20;
    double a = -2;
    double b = 2;
    double *x = new double[N];
    double *w = new double[N];
    // Set up mesh points and weights
    clock_t start, finish;

    start = clock();
    gauleg(a,b,x,w,N);
    // Evaluate the integral with Gauss-Legendre method

    double int_gauss = 0.;
    for (int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
            for (int k = 0; k<N; k++){
                for (int l = 0; l<N; l++){
                    for (int m = 0; m<N; m++){
                        for (int n = 0; n<N; n++){
        int_gauss += w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_function(x[i],x[j],x[k],x[l],x[m],x[n]);
                        }
                    }
                }
            }
        }
    }
    finish = clock( );
     ((finish - start)/CLOCKS_PER_SEC ) ;
    cout << "Finish" << endl;
    cout << finish << endl;
     cout << "start" << endl;
     cout << start << endl;
     cout << "Finish - Start" << endl;
     cout << finish - start << endl;
    cout << "Gaussian quad = " << int_gauss << endl;
    output(N, a, b, int_gauss);
    delete [] x;
    delete [] w;
    return 0;
}

 // Define the function to integrate: Six-dimension integral

double int_function(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double alpha = 2.;  // Alpha: Correspond to the charge of the helium atom Z=2
    // Evaluation of the different terms of the exponential
    double exp1 = -2*alpha*sqrt((x1*x1)+(y1*y1)+(z1*z1));
    double exp2 = -2*alpha*sqrt((x2*x2)+(y2*y2)+(z2*z2));
    double deno = sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
    if(deno < pow(10.,-6)) {return 0;}
    else return exp(exp1+exp2)/deno;
}
void output(int N, double a, double b, double int_gauss)
{

      ofstream writer( "Project3a.txt" , ios::app ) ;

      writer << "RESULTS Project 3a:" << endl;
      writer << setiosflags(ios::showpoint | ios::uppercase);
      writer <<"N = " << setw(15) << N << endl;
      writer <<"a = " << setw(15) << a << endl;
      writer <<"b = " << setw(15) << b << endl;
      writer <<"int_gauss = " << setw(15) << int_gauss << endl;
      writer.close() ;
}
