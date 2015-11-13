#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "../libraries/lib.h"
using namespace  std;
using namespace arma;

// output file as global variable
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit); // limit = storrelse grid // add = +-1
}

// Function declaration
void output(int n_spins, int mcs, double temperature, double *average);   // write to .txt function
void initialize(int n_spins, double temperature, int **spin_matrix, double& E, double& M);
void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w);


int main()
{
    // initial variables
    int mcs = 10000; // Number of Monte Carlo trials
    int n_spins = 80; // Lattice size or number of spins (x and y equal)
    int n_temp = 40; // Number of temperatures

    double initial_temp = 2.; //
    double final_temp = 2.4; //
    //double temp_step = 0.5; //
    double temp_step = (final_temp - initial_temp)/(n_temp-1);
    double temperature;

    int **spin_matrix;
    double w[17], average[5], E, M;

    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    long idum = -1; // random starting point

    // Lager vec for Matlab

    vec t = zeros<vec>(n_temp); // Temp
    vec AM = zeros<vec>(n_temp); //

    vec AE = zeros<vec>(n_temp);
    vec AE2 = zeros<vec>(n_temp);
    //vec AM = zeros<vec>(n_temp);
    vec AM2 = zeros<vec>(n_temp);
    vec AabsM = zeros<vec>(n_temp);

    vec Evariance = zeros<vec>(n_temp);
    vec Mvariance = zeros<vec>(n_temp);
    vec M2variance = zeros<vec>(n_temp);

    vec HeatCv = zeros<vec>(n_temp);
    vec Susceptibility = zeros<vec>(n_temp);

    int counter = 0;  // Used to count temp

    for ( temperature = initial_temp; temperature <= final_temp; temperature+=temp_step){
        //    initialise energy and magnetization
        E = M = 0.;
        // setup array for possible energy changes
        for( int de =-8; de <= 8; de++) w[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);

        // initialise array for expectation values
        for( int i = 0; i < 5; i++) average[i] = 0.;
        initialize(n_spins, temperature, spin_matrix, E, M);

        // start Monte Carlo computation
        for (int cycles = 1; cycles <= mcs; cycles++){
            Metropolis(n_spins, idum, spin_matrix, E, M, w);
            // update expectation values
            average[0] += E;    average[1] += E*E;
            average[2] += M;    average[3] += M*M; average[4] += fabs(M);
        }
        //cout << temperature << endl; // temperature range in loop
        t(counter) = temperature;
        AE(counter) = average[0]/mcs;
        AE2(counter) = average[1]/mcs;
        AM(counter) = average[2]/mcs;
        AM2(counter) = average[3]/mcs;
        AabsM(counter) = average[4]/mcs;

        Evariance(counter) = (AE2(counter)- AE(counter)*AE(counter))/n_spins/n_spins;
        Mvariance(counter) = (AM2(counter) - AM(counter)*AM(counter))/n_spins/n_spins;
        M2variance(counter) = (AM2(counter) - AabsM(counter)*AabsM(counter))/n_spins/n_spins;

        HeatCv(counter) = Evariance(counter)/t(counter)/t(counter);
        Susceptibility(counter) = M2variance(counter)/t(counter);

        counter ++; //

    }

    t.save("t.txt", raw_ascii);
    //AM.save("AM.txt", raw_ascii);
    AE.save("AE80.txt", raw_ascii);
    //AE2.save("AE2.txt", raw_ascii);
    //AM2.save("AM2.txt", raw_ascii);
    AabsM.save("AabsM80.txt", raw_ascii);

    //Evariance.save("Evariance.txt", raw_ascii);
    //Mvariance.save("Mvariance.txt", raw_ascii);
    //M2variance.save("M2variance.txt", raw_ascii);

    HeatCv.save("HeatCv80.txt", raw_ascii);
    Susceptibility.save("Susceptibility80.txt", raw_ascii);


    output(n_spins, mcs, temperature, average);     // Write to .txt
    cout << "Lattice size: " << n_spins << endl;
    cout << "Mcs: " << mcs << endl;
    cout << "Average: " << *average << endl;
    cout << "Temp " << temperature << endl;  // temperature after loop: temperature+=temp_step

    return 0;

}


void output(int n_spins, int mcs, double temperature, double *average)
{
    // .txt file start
    ofstream writer( "Prosjekt 4bV5.txt" , ios::app ) ;
    //   Read to file
    writer << "RESULTS Prosjekt 4b:" << endl;
    writer << "Lattice size: " << n_spins << endl;
    writer << "Mcs: " << mcs << endl;
//    writer << "Temperature: " << temperature << endl;
//    writer << "Average: " << *average << endl;
    writer << ""  << endl;

    double norm = 1/((double) (mcs)); // divided by total number of cycles
    double Eaverage = average[0]*norm;
    double E2average = average[1]*norm;
    double Maverage = average[2]*norm;
    double M2average = average[3]*norm;
    double Mabsaverage = average[4]*norm;

    // all expectation values are per spin, divide by 1/n_spins/n_spins
    double Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins;
    double Mvariance = (M2average - Maverage*Maverage)/n_spins/n_spins;
    double M2variance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;


    writer << "Temperature: " << temperature << endl;
    writer << "Average: " << *average << endl;
    writer << "Eaverage/n_spins/n_spins " << Eaverage/n_spins/n_spins << endl;
    writer << "HeatCv = Evariance/temperature/temperature " << Evariance/temperature/temperature << endl;
    writer << "Susceptibility = M2variance/temperature " << M2variance/temperature << endl;
    writer << "Mabsaverage/n_spins/n_spins " << Mabsaverage/n_spins/n_spins << endl;
    writer << ""  << endl;

    writer.close() ;

    // .txt file writer finished

}


// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, double temperature, int **spin_matrix,
                double& E, double& M)
{
//     setup spin matrix and intial magnetization
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
    // Step 1:Initial state: position random in lattice
    for(int i=0; i<n_spins*n_spins; i++) {
        int ix = (int) (ran1(&idum)*(double)n_spins);
        int iy = (int) (ran1(&idum)*(double)n_spins);

        int deltaE =  2*spin_matrix[iy][ix]*
                (spin_matrix[iy][periodic(ix,n_spins,-1)]+
                spin_matrix[periodic(iy,n_spins,-1)][ix] +
                spin_matrix[iy][periodic(ix,n_spins,1)] +
                spin_matrix[periodic(iy,n_spins,1)][ix]);  // Energy difference

        if ( ran1(&idum) <= w[deltaE+8] ) {     // Condition
            spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
            M += (double) 2*spin_matrix[iy][ix];
            E += (double) deltaE;
        }

    }
} // end of Metropolis sampling over spins

