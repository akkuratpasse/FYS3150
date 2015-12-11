#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "ran0.cpp"
#include "time.h"

using namespace std;
using namespace arma;
int nplanets = 50;
long idum=-1 ;
double e = 0.01;
double x,y,z;
vec F(vec V){

    vec F = zeros<vec>(3*nplanets);
    for (int i=0;i<nplanets;i++){
        F[3*i+0]=V[3*i+0];F[3*i+1]=V[3*i+1];F[3*i+2]=V[3*i+2];
    }
    return F;
}
vec E(vec X,vec M){

    vec E = zeros<vec>(3*nplanets);
    double grav=4*3.141592*3.141592;
    for (int j=0;j<nplanets;j++){
        for (int i=j+1;i<nplanets;i++){
            double dr=sqrt(pow((X[3*j+0]-X[3*i+0]),2) + pow((X[3*j+1]-X[3*i+1]),2) +
                    pow((X[3*j+2]-X[3*i+2]),2));
            double  f=-grav*M[i]*M[j]/((pow(dr,2)+e)*dr);
            E[3*i+0]+=(f/M[i])*(X[3*i+0]-X[3*j+0]);
            E[3*i+1]+=(f/M[i])*(X[3*i+1]-X[3*j+1]);
            E[3*i+2]+=(f/M[i])*(X[3*i+2]-X[3*j+2]);
            E[3*j+0]+=(f/M[j])*(X[3*j+0]-X[3*i+0]);
            E[3*j+1]+=(f/M[j])*(X[3*j+1]-X[3*i+1]);
            E[3*j+2]+=(f/M[j])*(X[3*j+2]-X[3*i+2]);
        }
    }
    return E;
}

vec K(vec V, vec M){
    vec K = zeros<vec>(nplanets);
    for(int i=0;i<nplanets;i++){
        K[i]=0.5*M[i]*M[i]*(V[3*i+0]*V[3*i+0]+V[3*i+1]*V[3*i+1]+V[3*i+2]*V[3*i+2]);
}

    return K;
}
vec P(vec X, vec M){
    vec P = zeros<vec>(nplanets);
    double grav=4*3.141592*3.141592;
    for (int j=0;j<nplanets;j++){
        for (int i=j+1;i<nplanets;i++){
            double dr=sqrt(pow((X[3*j+0]-X[3*i+0]),2) + pow((X[3*j+1]-X[3*i+1]),2) +
                    pow((X[3*j+2]-X[3*i+2]),2));
            double  f=-grav*M[i]*M[j]/dr;
        P[i]+=f;
        P[j]+=f;
    }
    return P;
}
}
int main()
{
    clock_t start,finish;
    int counter=0;
    int N=100 ;double T=0.1;
    vec X = zeros<vec>(3*nplanets);
    vec M = zeros<vec>(nplanets);
    vec V = zeros<vec>(3*nplanets);
    vec k = zeros<vec>(nplanets);
    vec p = zeros<vec>(nplanets);
    vec gold = zeros<vec>(3*nplanets);
    vec g = zeros<vec>(3*nplanets);
    mat Z = zeros<mat>(3*nplanets,N+1);
    double t=0.0,dt=T/N;
//    M[0]=1;M[1]=pow(10,-8);  // For the 2-body problem, also comment the lines 94-97
//    X[0]=0.0;X[1]=0.0;X[2]=0.0;X[3]=0.0;X[4]=1.0;X[5]=0.0;
//    V[0]=0.0;V[1]=0.0;V[2]=0;V[3]=2*3.14;V[4]=0.0;V[5]=0.0;
    for (int i=0;i<nplanets;i++){
        int ro=20;
        bool not_inside=true;
        while(not_inside){
             x=2*ro*ran0(&idum)-ro;
             y=2*ro*ran0(&idum)-ro;
             z=2*ro*ran0(&idum)-ro;
            if((x*x+y*y+z*z)<ro*ro){
                not_inside=false;
            }
        }
        X(3*i+0)=x;
        X(3*i+1)=y;
        X(3*i+2)=z;
        M(i)=1;
    }
    gold=E(X,M);
    k = K(V,M);
    p = P(X,M);
    cout << "Total Energy (before)" << sum(k+(p/2)) << endl;
    start = clock();
    while(t <= T){
        X = X + V*dt+0.5*gold*dt*dt;
        g=E(X,M);
        V = V + 0.5*(g+gold)*dt;
        gold=g;
        t = t+dt;
        Z.col(counter)=X;
        counter++;
        //cout << K << endl;
    }
    finish = clock();
    ((finish-start)/ CLOCKS_PER_SEC);
    k = K(V,M);
    p = P(X,M);
    Z.save("Z.txt", raw_ascii);
      cout << "Total Time" << finish-start << endl;
      cout << "Start Time" << start << endl;
      cout << "Finish Time" << finish << endl;
      //cout << "Kinetic Energy" << k << endl;
     //cout << "Potential Energy" << p << endl;
      cout << "Total Energy (after)" << sum(k+(p/2)) << endl;
      cout << 2*sum(k)<< endl;
      cout << -sum(p/2) << endl;
  //  cout << "X" << K << endl;
   // cout << "F" << V << endl;
//    cout << "t" << t << endl;
    return 0;
}

