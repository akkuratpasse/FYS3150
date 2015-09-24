#include <iostream>
#include <armadillo>
#include "../libraries/lib.h"
using namespace  std;
using namespace arma;

// Function declaration
double offdiag(mat A, int &p, int &q, int n);
void JacobiRotate(mat &A, mat &R, int k, int l, int n);


// Define a matrix A and a matrix R for the eigenvector

int n = 4;
int j = n;
int i = n;

int main(int argc, char** argv){
    // Creating matrix A and identity matrix R
    //mat X;  //  X test
    //X << 2 << 1 << endr
    //  << 2 << 1 << endr;

    mat I = eye<mat>(i,j); // Identity matrix with Armadillo
    cout << "I:" << endl << I << endl;

    mat A = zeros<mat>(n,n);
    A.diag() += 2;      // pa diagonalen
    A.diag(-1) += 1;    // rett under diagonalen
    A.diag(1) += 1;     // rett over diagonalen
    cout << "A:" << endl << A << endl;

    mat R = randu<mat>(i,j); // Custom Identity matrix R
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            if (i == j) {
                R(i,j)=1.;
               } else{
                R(i,j)=0;
                }
            }
    }
    cout << "R:" << endl << R << endl;



    //
    double tolerance = 1.0E-10;
    int iterations = 0;  // total iterations
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

    cout << "A:" << endl << A << endl;

    return 0;
}

// Offdiagonal function
double offdiag(mat A, int &p, int &q, int n)
{
    double max;
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
        tau = (A(1,1) - A(k,k))/(2*A(k,1));

        if (tau >= 0){
            t = 1.0/(tau +sqrt(1.0 + tau*tau));
        } else {
            t = -1.0/(-tau +sqrt(1.0 + tau*tau));
        }
        c = 1/sqrt(1+t*t);
        s = c*t;
    } else {
        c = 1.0;
        s = 1.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;   //
    A(l,k) = 0.0;   //
    for (int i = 0; i < n; i++){
        if (i != k && i != 1){
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
