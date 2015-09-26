/*
  Solves the one-particle Schrodinger equation
*/

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

double potential(double);
int comp(const double *, const double *);
void output(double r_min , double r_max, int max_step, double *d);

// initialise constants
double r_min = 0;
double r_max = 5;
int orb_l = 0;
int max_step = 10;
double step    = (r_max - r_min) / max_step;
double const_2 = -1.0 / (step * step);
double const_1 =  - 2.0 * const_2;
double orb_factor = orb_l * (orb_l + 1);


// Start Main
int main(int argc, char* argv[])
{
  int       i, j;
  double    *e, *d, *w, *r, **z;

  // local memory for r and the potential w[r]
  r = new double[max_step + 1];
  w = new double[max_step + 1];
  for(i = 0; i <= max_step; i++) {
    r[i] = r_min + i * step;
    w[i] = potential(r[i]) + orb_factor / (r[i] * r[i]);
  }
  // local memory for the diagonalization process
  d = new double[max_step];    // diagonal elements
  e = new double[max_step];    // tri-diagonal off-diagonal elements
  z = (double **) matrix(max_step, max_step, sizeof(double));
  for(i = 0; i < max_step; i++) {
    d[i]    = const_1 + w[i + 1];
    e[i]    = const_2;
    z[i][i] = 1.0;
    for(j = i + 1; j < max_step; j++)  {
      z[i][j] = 0.0;
    }
  }
  // diagonalize and obtain eigenvalues
  tqli(d, e, max_step - 1, z);
  // Sort eigenvalues as an ascending series
  qsort(d,(UL) max_step - 1,sizeof(double),
         (int(*)(const void *,const void *))comp);


    // Write to file function
    output(r_min , r_max, max_step, d);

    return 0;
} // End: function main()

// End Main

/*
  The function potential()
  calculates and return the value of the
  potential for a given argument x.
  The potential here is for the hydrogen atom
*/

double potential(double x)
{
   return -2./x;

} // End: function potential()


/*
  The function   int comp()
  is a utility function for the library function qsort()
  to sort double numbers after increasing values.
*/

int comp(const double *val_1, const double *val_2)
{
  if((*val_1) <= (*val_2))       return -1;
  else  if((*val_1) > (*val_2))  return +1;
  else                     return  0;
} // End: function comp()


void output(double r_min , double r_max, int max_step, double *d)
{
    // .txt file start
      ofstream writer( "Householder.txt" , ios::app ) ;

//      if( ! writer )
//      {
//        cout << "Error opening file for output" << endl ;
//        return -1 ;
//      }

      //   Read to file

  //    writer << "RESULTS Project 2 Householder Algorithm:" << endl;
  //    writer << setiosflags(ios::showpoint | ios::uppercase);
  //    writer <<"N = " << setw(15) << n << endl;
  //    writer <<"R min = " << setw(15) << setprecision(4) << r_min << endl;
  //    writer <<"R max = " << setw(15) << setprecision(4) << r_max << endl;
  //    writer << "Number of iterations:" << setw(5) << iterations << endl;
  //    writer << "Three lowest eigenvalues from rotated matrix:" << endl;
  //    writer << "B[0]  = " << setw(15) << setprecision(5) << B[0] << endl;
  //    writer << "B[1]  = " << setw(15) << setprecision(5) << B[1] << endl;
  //    writer << "B[2]  = " << setw(15) << setprecision(5) << B[2] << endl;
  //    writer << "" << endl;


      writer << "RESULTS Project 2 Householder Algorithm:" << endl;
      writer << setiosflags(ios::showpoint | ios::uppercase);
      writer <<"R_min = " << setw(15) << setprecision(8) << r_min << endl;
      writer <<"R_max = " << setw(15) << setprecision(8) << r_max << endl;
      writer <<"Number of steps = " << setw(15) << max_step << endl;
      writer << "Five lowest eigenvalues:" << endl;
      for(int i = 0; i < 5; i++) {
        writer << setw(15) << setprecision(8) << d[i] << endl;
      }

      writer.close() ;

    // .txt file writer finished
}
