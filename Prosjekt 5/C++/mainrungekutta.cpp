#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "ran0.cpp"
#include "time.h"

using namespace std;
using namespace arma;
int nplanets = 100;
long idum=-1 ;
double e = 0;
double x,y,z;
//vec F(vec X,vec M){

//    vec F = zeros<vec>(6*nplanets);
//    double grav=4*3.141592*3.141592, msun=1;//datum::pi
//    for (int i=0;i<nplanets;i++){
//    F[6*i+0]=X[4*i+2];F[4*i+1]=X[4*i+3];
//    F[4*i+2]=-grav*msun*X[4*i+0]/pow(((X[4*i+0]*X[0]) + (X[4*i+1]*X[4*i+1])),3.0/2);
//    F[4*i+3]=-grav*msun*X[4*i+1]/pow(((X[4*i+0]*X[0]) + (X[4*i+1]*X[4*i+1])),3.0/2);

//    }
//    return F;
//}
vec E(vec X,vec M){

    vec E = zeros<vec>(6*nplanets);
    double grav=4*3.141592*3.141592;
    for (int j=0;j<nplanets;j++){
        E[6*j+0]=X[6*j+3];E[6*j+1]=X[6*j+4];E[6*j+2]=X[6*j+5];
        for (int i=j+1;i<nplanets;i++){
            double dr=sqrt(pow((X[6*j+0]-X[6*i+0]),2) + pow((X[6*j+1]-X[6*i+1]),2) +
                     pow((X[6*j+2]-X[6*i+2]),2));
            //dr=sqrt(pow((X[3*j+0]-X[3*i+0]),2) + pow((X[3*j+1]-X[3*i+1]),2) +
              //                  pow((X[3*j+2]-X[3*i+2]),2));
            double  f=-grav*M[i]*M[j]/((pow(dr,2)+e)*dr);
            E[6*i+3]+=(f/M[i])*(X[6*i+0]-X[6*j+0]);
            E[6*i+4]+=(f/M[i])*(X[6*i+1]-X[6*j+1]);
            E[6*i+5]+=(f/M[i])*(X[6*i+2]-X[6*j+2]);
            E[6*j+3]+=(f/M[j])*(X[6*j+0]-X[6*i+0]);
            E[6*j+4]+=(f/M[j])*(X[6*j+1]-X[6*i+1]);
            E[6*j+5]+=(f/M[j])*(X[6*j+2]-X[6*i+2]);
        }
    }

    //cout << X << endl;
   // exit(1);
    return E;
}
vec K(vec X, vec M){
    vec K = zeros<vec>(nplanets);
    for(int i=0;i<nplanets;i++){
        K[i]=0.5*M[i]*(X[6*i+3]*X[6*i+3]+X[6*i+4]*X[6*i+4]+X[6*i+5]*X[6*i+5]);
    }
    return K;
}
vec P(vec X, vec M){
    vec P = zeros<vec>(nplanets);
    double grav=4*3.141592*3.141592;
    for (int j=0;j<nplanets;j++){
        for (int i=j+1;i<nplanets;i++){
            double dr=sqrt(pow((X[6*j+0]-X[6*i+0]),2) + pow((X[6*j+1]-X[6*i+1]),2) +
                    pow((X[6*j+2]-X[6*i+2]),2));
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
    int N=100; double T=10;
    vec X = zeros<vec>(6*nplanets);
    vec M = zeros<vec>(nplanets);
    vec V = zeros<vec>(6*nplanets);
    mat Z = zeros<mat>(6*nplanets,N+1);
    vec k = zeros<vec>(nplanets);
    vec p = zeros<vec>(nplanets);
    double t=0.0,dt=T/N;
    vec F1 = zeros<vec>(6*nplanets);
    vec F2 = zeros<vec>(6*nplanets);
    vec F3 = zeros<vec>(6*nplanets);
    vec F4 = zeros<vec>(6*nplanets);
//    M[0]=1;M[1]=pow(10,-8);
//    X[0]=1.0;X[1]=0.0;X[2]=0.0;X[3]=0.0;X[4]=0.0;X[5]=0.0;
//    V[0]=2*3.14;V[1]=0.0;V[2]=0;V[3]=0.0;V[4]=0.0;V[5]=0.0;
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
        X(6*i+0)=x;
        X(6*i+1)=y;
        X(6*i+2)=z;
        M(i)=1;
    }
    k=K(X,M);
    p=P(X,M);
    F1=E(X,M);
   // cout << k << endl;
   // cout << p << endl;
   // cout << sum(k+(p/2)) << endl;
    start = clock();
    while(t <= T){

        X = X + V*dt;
        F1=E(X,M); F2=E(X+(dt/2)*F1,M);
        F3=E(X +(dt/2)*F2,M);
        F4=E(X + dt*F2,M);
        X=X + (dt/6)*(F1+2*F2+2*F3+F4);
        if (isnan(X(1))){
            cout << "n a n"<< endl;
        }
        t = t+dt;
        Z.col(counter)=X;
        counter++;

        //cout << X << endl;


}
    k=K(X,M);
    //cout << sum(k+(p/2)) << endl;
   // cout << k << endl;
    finish = clock();
//    for(int j=0;j<1;j++){
//        Z[j]=x;Z[j+1]=y;Z[j+2]=vx;Z[j+3]=vy;
//    }
    ((finish-start)/ CLOCKS_PER_SEC);
    Z.save("Z.txt", raw_ascii);
      cout << finish-start << endl;
      cout << start << endl;
      cout << finish << endl;
//    cout << "X" << X << endl;
//    cout << "F" << V << endl;
//    cout << "t" << t << endl;
    return 0;
}

