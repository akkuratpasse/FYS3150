#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <iomanip>
#include "../libraries/lib.h"
using namespace  std;
using namespace arma;
// output file as global variable
ofstream ofile;

// Function declaration
double offdiag(mat A, int &p, int &q, int n);
void JacobiRotate(mat &A, mat &R, int k, int l, int n);
void output(int n, double rmin, double rmax, double wr, int iterations, mat &B);   // write to .txt function

// Define a matrix A and a matrix R for the eigenvector

int n = 100; // Size of matrix    // Need to change
int j = n;
int i = n;

int main(int argc, char** argv){

 // Equation
    double rmax = 10.;   //Max values  // Need to change
    double rmin = 0;   //Min values
    double h = (rmax-rmin)/(n+1);   // Step length
    double wr = 0.01;   // Parameter strength of the oscillator potential

    // Define Diagonal matrix element, d
    mat A = zeros<mat>(n,n);
    vec r = zeros<vec>(n);
    vec v = zeros<vec>(n);
    for (int i = 0; i<n; i++){
        r[i]= rmin + (i+1)*h;
        //v[i]= r[i]*r[i]; // Potential
        v[i] = wr*wr*r[i]*r[i]+(1/r[i]);    // New Potential
        A.diag()[i] += 2/(h*h)+v[i];     // on the diagonal
    }
    //cout << "ro, Dimensionless variable" << endl;
    //cout << r << endl;
    //cout << "v, Harmonic oscillator potential" << endl;
    //cout << v << endl;

    // Define Non-diagonal matrix element, e
    double e = -(1/(h*h));
    A.diag(-1) += e;    // rett below diagonalen
    A.diag(1) += e;    // rett over diagonalen

    cout << "A:" << endl << A << endl;
    //A.save("A.txt", raw_ascii);

    // Armadillo function finding Eigenvalue for mat A
//    vec E = zeros<vec>(n);  // lage vektor for Eigenvalue
//    eig_sym( E, A );

//    cout << "E0 " << E[0] << endl;      //E0 1.24981 n=100
//    cout << "E1 " << E[1] << endl;      //E1 2.18925 n=100

    // Armadillo function finding Eigenvalue, Eigenvector for mat A)
    vec E = zeros<vec>(n);  // Armadillo, lage vektor for Eigenvalue
    mat Z = zeros<mat>(n,n);  // Armadillo, lage mat for Eigenvektor
    eig_sym(E, Z, A);
    cout << "Armadillo E0 " << E[0] << endl;      //E0 1.24981 n=100
    cout << "Armadillo E1 " << E[1] << endl;      //E1 2.18925 n=100
    cout << "Armadillo E2 " << E[2] << endl;
    cout << "Armadillo Eigenvektor V0 " << Z.col(0) << endl;
//    cout << "V1 " << Z.col(1) << endl;
//    E.save("Earmadillo.txt", raw_ascii);
//    Z.save("Zarmadillo.txt", raw_ascii);

    // Save files for Plot
    Z.save("ZArmadilloN100_Omega001.txt", raw_ascii);    // Armadillo Eigenvectors used to plot
    r.save("rhoN100_Omega001.txt", raw_ascii);           // rho used to plot

    mat R = zeros<mat>(i,j); // Custom Identity matrix R
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            if (i == j) {
                R(i,j)=1.;
               } else{
                R(i,j)=0;
                }
            }
    }

    double tolerance = 1.0E-10;
    int iterations = 0;  // Iterations
    double maxnondiag = offdiag(A, i, j ,n);
    double maxiter = (double) n * (double) n * double (n); // 64
    while (maxnondiag > tolerance && iterations <= maxiter)
    {
        int p, q;
        // Use functions
        maxnondiag = offdiag(A, i, j ,n);
        JacobiRotate(A,R,i,j,n);
        iterations++;
    }
    std:: cout << "Number of iterations:" << iterations << "/n" << endl;

    cout << "ARotated:" << endl << A << endl;
    //A.save("ARotated.txt", raw_ascii);

    // Eigenvalues
    mat B = sort(A.diag());     // Eigenvalues Sorted Low to high
    cout << "Eigenvalues, sorted low to high after " << iterations << " Jacobi Rotations" << endl;
    for(int i = 0; i < 3; i++) {
      cout << setw(15) << setprecision(6) << B[i] << endl;
    }


    output(n,rmin, rmax, wr, iterations, B);     // Write to .txt

}

// Offdiagonal function
double offdiag(mat A, int &p, int &q, int n)
{
    double max=0;
    for (int i=0; i<n; i++){
        for (int j=i+1; j<n; j++){
        {
            double aij = fabs(A(i,j));
            if (aij > max)
            {
                max = aij;
                p = i;
                q = j;
                }
            }
        }
    }
    //cout << "From function offdiag() Max value is:" << endl << max << endl;
    return max;
}

// Jacobi_rotate
void JacobiRotate(mat &A, mat &R, int k, int l, int n)
{
    double s, c;
    if (A(k,l) != 0.0){
        double t, tau;
        tau = (A(l,l) - A(k,k))/(2*A(k,l)); // Tau

        if (tau > 0){
            t = 1.0/(tau +sqrt(1.0 + tau*tau));
        } else {
            t = -1.0/(-tau +sqrt(1.0 + tau*tau));
        }
        c = 1/sqrt(1+t*t);
        s = c*t;
    } else {
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;

    double A_kl = (a_kk - a_ll)*c*s + A(k,l)*(c*c-s*s);
    A(k,l) = 0.0;   //
    A(l,k) = 0.0;   //
    for (int i = 0; i < n; i++){
        if (i != k && i != l){
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
    // New Eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;

    }
    return;
}

// Write to file
void output(int n, double rmin, double rmax, double wr, int iterations, mat &B)
{
    // .txt file start
      ofstream writer( "HarmonicOscillatorPotential100_Omega001.txt" , ios::app ) ;

      //   Read to file

      writer << "RESULTS Project 2c Harmonic Oscillator Potential" << endl;
      writer << setiosflags(ios::showpoint | ios::uppercase);
      writer << "N = " << setw(15) << n << endl;
      writer << "R min = " << setw(15) << setprecision(4) << rmin << endl;
      writer << "R max = " << setw(15) << setprecision(4) << rmax << endl;
      writer << "Strenght of the Oscillator Potential = " << setw(5) << setprecision(2) << wr << endl;
      writer << "Number of iterations:" << setw(5) << iterations << endl;
      writer << "Three lowest eigenvalues from rotated matrix:" << endl;
      for(int i = 0; i < 3; i++) {
        writer << setw(15) << setprecision(8) << B[i] << endl;
      }
      writer << "" << endl;

      writer.close() ;

    // .txt file writer finished
}
