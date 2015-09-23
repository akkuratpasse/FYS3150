#include <iostream>
#include <armadillo>
#include "../libraries/lib.h"
using namespace  std;
using namespace arma;

// Function declaration
double offdiag(mat A, int p, int q, int n);

// Define a matrix A and a matrix R for the eigenvector

int n = 2;
int j = n;
int i = n;

int main(int argc, char** argv){
    // Creating matrix A and identity matrix R
    mat X;
    X << 2 << 1 << endr
      << 2 << 1 << endr;
   // mat A = randu<mat>(i,j);
    mat R = eye<mat>(i,j);

    mat A = randu<mat>(i,j);
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            if (i == j) {
                A(i,j)=1.;
               } else{
                A(i,j)=0;
                }
            }
}
    cout << "X:" << endl << X << endl;
    cout << "R:" << endl << R << endl;
    cout << "A:" << endl << A << endl;

    // Use functions
    offdiag(A, n, n ,n);


    return 0;
}

// Offdiagonal function
double offdiag(mat A, int p, int q, int n)
{
    double max;
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
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
    cout << "From function offdiag() Max value is:" << endl << max << endl;
    return max;
}


