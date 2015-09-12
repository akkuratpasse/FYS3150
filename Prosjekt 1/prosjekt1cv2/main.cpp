#include <math.h>
#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;

// Project 1 - Gaussian Elimination and Tridiagonal matrices
// Part 2 - Relative Error med forskjellig h

vec vsolver(int n)  // Funksjon - tridiagonal system solver
{
    // Numerical Solution
    double h = 1.0/ ((double) n+1);      // Step length

    // Use Armadillo to allocate arrays
    vec a = (-1)*ones<vec>(n);
    vec b = 2*ones<vec>(n);
    vec c = (-1)*ones<vec>(n);

    vec v = zeros<vec>(n);
    vec f = zeros<vec>(n);

//    cout << a << endl;        // Hva som skal streames
//    cout << b << endl;
//    cout << c << endl;


    // Define value of righthand side
    // b = h*h*fi, fi = 100*exp^(-10*x)
    for(int i=0; i < n; i++){
        f[i] = h*h*(100*exp(-10*h*(i+1)));
    }
//    cout << "f, righthand side" << endl;
//    cout << f << endl;


    // Define value of lefthandside

    // Forward Substitution
    for(int i=1; i < n; i++){
        double m = a[i-1]/b[i-1];
        b[i] = b[i]-(m*c[i-1]);
        f[i] = f[i]-(m*f[i-1]);
    }
//    cout << "Forward Substitution" << endl;
//    cout << "b" << endl;
//    cout << b << endl;
//    cout << "f" << endl;
//    cout << f << endl;


    // hva siste b er i diagonalen
    v[n-1] = f[n-1]/b[n-1];
    //cout << v[n-1] << endl; // Skriver ut siste tall pa vektor

    // Backward Substitution
    for(int i = n-1; i > 0; i--){
        v[i-1] = (f[i-1]- c[i-1]*v[i])/b[i-1];
    }
//      cout << "Backward Substitution" << endl;
//    cout << "v, Numerical Solution" << endl;
//    cout << v << endl;  // det som skrives ut er internal points

    v.save("V.txt", raw_ascii);

    // Analytical Solution
    // int n = 10;
    // Step length
    //  double h = 1.0/ ((double) n+1);

    // Use Armadillo to allocate arrays
    vec u = zeros<vec>(n);

    // Analytical Solution, u
    for(int i=0; i < n; i++){
        u[i] = (1-(1-exp(-10))*((i+1)*h)-exp(-10*(i+1)*h));
    }
//    cout << "u, Analytical Solution" << endl;
//    cout << u << endl;  // det som skrives ut er internal points
//    u.save("U.txt", raw_ascii);

    return v;
}

int main(int argc, char** argv)
{
    // Relative Error med forskjellig h

    int n, j;
    double h;
    // Deklarer vektor
    vec v = zeros<vec>(n);
    vec e2 = zeros<vec>(6);
    for(int i=0; i < 6; i++){
        n = int (pow(10,i)); //
        h = 1.0/(n+1);   // forskjellige h er avhengig av n
        v = vsolver(n);     // Henter numerisk funksjon
        vec u = zeros<vec>(n);      // Vektor for Analytisk
        for(int k=0; k < n; k++){
            u[k] = (1-(1-exp(-10))*((k+1)*h)-exp(-10*(k+1)*h)); // Analytisk
        }
//        cout << "u, Analytisk" << endl;
//        cout << u << endl;
        j=n/2;      // Velge ett fast punkt
        e2[i] = log10(abs(((v[j]-u[j])/u[j])));     // Relative Error
    }
    cout << "e2, Relative Error log10 with different h" << endl;
    cout << e2 << endl;
    e2.save("e2.txt");

    return 0;
}
