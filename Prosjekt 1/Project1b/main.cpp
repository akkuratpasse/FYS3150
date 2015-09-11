#include <math.h>
#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;

// Project 1 - Gaussian Elimination and Tridiagonal matrices

int main(int argc, char** argv)
{
    int n = 10;

    // Step length
    double h = 1.0/ ((double) n+1);

    // Use Armadillo to allocate arrays
    vec a = (-1)*ones<vec>(n);
    vec b = 2*ones<vec>(n);
    vec c = (-1)*ones<vec>(n);

    vec v = zeros<vec>(n);
    vec f = zeros<vec>(n);

    cout << a << endl;
    cout << b << endl;
    cout << c << endl;


    // Define value of righthand side   
    // b = h*h*fi, fi = 100*exp^(-10*x)
    for(int i=0; i < n; i++){
        f[i] = h*h*(100*exp(-10*h*(i+1)));
    }
    cout << f << endl;


    // Define value of lefthandside
    // Ax = a(i)(i-1)*v(i-1)+ a(i)(i+1)*v(i+1)+a(i)(i)*v(i)

    // Forward Substitution
    for(int i=1; i < n; i++){
        double m = a[i-1]/b[i-1];
        b[i] = b[i]-(m*c[i-1]); // b starter pa b2 https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
        f[i] = f[i]-(m*f[i-1]);
    }
    cout << "b" << endl;
    cout << b << endl;
    cout << "f" << endl;
    cout << f << endl;


    // siden a diagonal na er 0 er det mulig a finne hva siste b er i diagonalen
    v[n-1] = f[n-1]/b[n-1];
//    cout << v[n-1] << endl; // Skriver ut siste tall pa diagonalen
//    v[n-1] = f[n-1]/b[n-1]; // versjon 2

    // Backward Substitution
    for(int i = n-1; i > 0; i--){
        v[i-1] = (f[i-1]- c[i-1]*v[i])/b[i-1];
    }

    cout << "v" << endl;
    cout << v << endl;
    // det som skrives ut er internal points
    v.save("V.txt", raw_ascii);

    //
}
