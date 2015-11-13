//#include <mpi.h>
//#include <iostream>
//using namespace std;

//int main(int args, char *argv[]) {
//    int numprocs;
//    int my_rank;
//    MPI_Init (&args, &argv);
//    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
//    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
//    cout << "I am processor " << my_rank << " of " << numprocs << endl;
//    MPI_Finalize();
//    return 0;
//}


/*
   Program to solve the two-dimensional Ising model
   with zero external field using MPI
   The coupling constant J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis sampling is used. Periodic boundary conditions.
   The code needs an output file on the command line and the variables mcs, nspins,
   initial temp, final temp and temp step.
*/
#include <mpi.h>
//#include "../libraries/mpi.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
using namespace  std;

// output file
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void initialize(int, int **, double&, double&);
// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *);
// prints to file the results of the calculations
void output(int, int, double, double *);
//  Matrix memory allocation
//  allocate space for a matrix
void  **matrix(int, int, int);
//  free space for  a matrix
void free_matrix(void **);
// ran2 for uniform deviates, initialize with negative seed.
double ran2(long *);

// Main program begins here

int main(int argc, char* argv[])
{
  char *outfilename;
  long idum;
  int **spin_matrix, n_spins, mcs, my_rank, numprocs;
  double w[17], average[5], total_average[5],
         initial_temp, final_temp, E, M, temp_step;

  //  MPI initializations
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0 && argc <= 1) {
    cout << "Bad Usage: " << argv[0] <<
      " read output file" << endl;
    exit(1);
  }
  if (my_rank == 0 && argc > 1) {
    outfilename=argv[1];
    ofile.open(outfilename);
  }
  n_spins = 2; mcs = 10000;  initial_temp = 1.; final_temp = 1; temp_step =0.5;
  /*
  Determine number of intervall which are used by all processes
  myloop_begin gives the starting point on process my_rank
  myloop_end gives the end point for summation on process my_rank
  */
  int no_intervalls = mcs/numprocs;
  int myloop_begin = my_rank*no_intervalls + 1;
  int myloop_end = (my_rank+1)*no_intervalls;
  if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = mcs;

  // broadcast to all nodes common variables
  MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //  Allocate memory for spin matrix
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  // every node has its own seed for the random numbers, this is important else
  // if one starts with the same seed, one ends with the same random numbers
  idum = -1-my_rank;  // random starting point
  // Start Monte Carlo sampling by looping over T first
  for ( double temperature = initial_temp; temperature <= final_temp; temperature+=temp_step){
    cout << "Started computing temperature: " << temperature << endl;
      //    initialise energy and magnetization
    E = M = 0.;
    // initialise array for expectation values
    initialize(n_spins, spin_matrix, E, M);
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
    for( int i = 0; i < 5; i++) average[i] = 0.;
    for( int i = 0; i < 5; i++) total_average[i] = 0.;
    // start Monte Carlo computation
    for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
      Metropolis(n_spins, idum, spin_matrix, E, M, w);
      // update expectation values  for local node
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);
    }

    // Find total average
    for( int i =0; i < 5; i++){
      MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    // print results
    if ( my_rank == 0) {
      output(n_spins, mcs, temperature, total_average);
    }
  }
  free_matrix((void **) spin_matrix); // free memory
  ofile.close();  // close output file
  // End MPI
  MPI_Finalize ();
  //return 0;
  cout << average[0] << endl;
}


// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, int **spin_matrix,
        double& E, double& M)
{
  // setup spin matrix and intial magnetization
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      spin_matrix[y][x] = 1; // spin orientation for the ground state
      M +=  (double) spin_matrix[y][x];
    }
  }
  // setup initial energy
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      E -=  (double) spin_matrix[y][x]*
    (spin_matrix[periodic(y,n_spins,-1)][x] +
     spin_matrix[y][periodic(x,n_spins,-1)]);
    }
  }
}// end function initialise

void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w)
{
  // loop over all spins
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      int ix = (int) (ran2(&idum)*(double)n_spins);
      int iy = (int) (ran2(&idum)*(double)n_spins);
      int deltaE =  2*spin_matrix[iy][ix]*
    (spin_matrix[iy][periodic(ix,n_spins,-1)]+
     spin_matrix[periodic(iy,n_spins,-1)][ix] +
     spin_matrix[iy][periodic(ix,n_spins,1)] +
     spin_matrix[periodic(iy,n_spins,1)][ix]);
      if ( ran2(&idum) <= w[deltaE+8] ) {
    spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
        M += (double) 2*spin_matrix[iy][ix];
        E += (double) deltaE;
      }
    }
  }
} // end of Metropolis sampling over spins


void output(int n_spins, int mcs, double temperature, double *total_average)
{
  double norm = 1/((double) (mcs));  // divided by total number of cycles
  double Etotal_average = total_average[0]*norm;
  double E2total_average = total_average[1]*norm;
  double Mtotal_average = total_average[2]*norm;
  double M2total_average = total_average[3]*norm;
  double Mabstotal_average = total_average[4]*norm;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2total_average- Etotal_average*Etotal_average)/n_spins/n_spins;
  double Mvariance = (M2total_average - Mabstotal_average*Mabstotal_average)/n_spins/n_spins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << Etotal_average/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
  ofile << setw(15) << setprecision(8) << Mtotal_average/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Mvariance/temperature;
  ofile << setw(15) << setprecision(8) << Mabstotal_average/n_spins/n_spins << endl;
} // end output function



//output   n_spins = 80; mcs = 1000000;  initial_temp = 2.1; final_temp = 2.6; temp_step =0.05;
//2.1000000     -1.6623223     0.96549405     0.86890125     0.94755969     0.86890125
//2.1500000     -1.6096958      1.1431762     0.83559545      1.8476218     0.83559545
//2.2000000     -1.5466559      1.4096649     0.78535838      4.4735811     0.78535838
//2.2500000     -1.4632705      2.0033194     0.64812222      31.390048     0.67402142
//2.3000000     -1.3504952      2.1416354    0.016253008      102.08798     0.37693760
//2.3500000     -1.2660605      1.4096325  -0.0026297728      48.923651     0.19094320
//2.4000000     -1.2037806      1.1360196   0.0049764538      22.320335     0.12434044
//2.4500000     -1.1520729     0.99399796  0.00080812594      13.912567    0.097844894
//2.5000000     -1.1062226     0.88016583 -5.8715625E-06      9.3244574    0.079459113
//2.5500000     -1.0651208     0.78682428 -0.00020778500      6.7595405    0.068241839
//2.6000000     -1.0283027     0.71879826 -0.00051571188      5.1270846    0.059789157

