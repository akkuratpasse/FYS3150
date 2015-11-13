#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cstdlib>
#include "../libraries/lib.h"
//#include "ran0.cpp"
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
void initialize(long &idum, int n_spins, double temperature, imat& spin_matrix, double& E, double& M);
//void initializeOrderDisorder(int spinorientation, int n_spins, long& idum, double temperature, imat &spin_matrix, double& E, double& M);
void Metropolis(int n_spins, long& idum, imat& spin_matrix, double& E, double&M, double *w);

int main()

{
    // initial variables

    double mcs = 10000; // Number of Monte Carlo trials
    int n_spins = 20; // Lattice size or number of spins (x and y equal)
    int n_temp = 1; // Number of temperatures
//    int spinorientation = 0; // 1=order, 0=dissorder

    //Temperatures
//    double initial_temp = 1.6; //
//    double final_temp = 2.6; //
    // low temp
//    double initial_temp = 1; //
//    double final_temp = 2.6;
    // high temp
    double initial_temp = 2.4;
    double final_temp = 2.4;

    //double temp_step = 0.1; //
    double temp_step = (final_temp - initial_temp)/(n_temp-1);
    double temperature;

    imat spin_matrix;
    double w[17], average[5], E, M;

    spin_matrix = zeros<imat>(n_spins,n_spins);
    long idum = -1; // random starting point

    // Lager vec for Matlab

    vec t = zeros<vec>(n_temp); // Temp
    vec AM = zeros<vec>(n_temp); //
//    vec AE = zeros<vec>(n_temp);
    vec AE2 = zeros<vec>(n_temp);
    //vec AM = zeros<vec>(n_temp);
    vec AM2 = zeros<vec>(n_temp);
    vec AabsM = zeros<vec>(n_temp);
    vec Evariance = zeros<vec>(n_temp);
    vec Mvariance = zeros<vec>(n_temp);
    vec M2variance = zeros<vec>(n_temp);
    vec HeatCv = zeros<vec>(n_temp);
    vec Susceptibility = zeros<vec>(n_temp);
//corder
    int tmp = (int) round(mcs);
    vec Mcs = zeros<vec>(tmp);
    vec AE = zeros<vec>(tmp);

    int counter = 0;  // Used to count temp

    for ( temperature = initial_temp; temperature <= final_temp; temperature+=temp_step){
        //    initialise energy and magnetization
        E = M = 0.;
        // setup array for possible energy changes
        for( int de =-8; de <= 8; de++) w[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);

        // initialise array for expectation values
        for( int i = 0; i < 5; i++) average[i] = 0.;
        initialize(idum, n_spins, temperature, spin_matrix, E, M);

        vec energies = zeros<vec>(mcs);
        vec magnetizations = zeros<vec>(mcs);
        // start Monte Carlo computation
        for (int cycles = 1; cycles <= mcs; cycles++){
            Metropolis(n_spins, idum, spin_matrix, E, M, w);
            // update expectation values
            average[0] += E;    average[1] += E*E;
            average[2] += M;    average[3] += M*M; average[4] += fabs(M);
            energies(cycles-1) = E;
            magnetizations(cycles-1) = M;
            Mcs(cycles-1) = cycles;
            AE(cycles-1) = average[0];
        }
        char filename[10000];
        sprintf(filename, "DISORDERenergiesT=%f.txt", temperature);
        energies.save(filename, raw_ascii);
        //cout << "Mcs" << Mcs << endl;
//        char filename[10000];
//        sprintf(filename, "OrdermagnetizationsT=%f.txt", temperature);
//        magnetizations.save(filename, raw_ascii);

        //cout << temperature << endl; // temperature range in loop
        t(counter) = temperature;
//        AE(counter) = average[0]/mcs;

        AE2(counter) = average[1]/mcs;
        AM(counter) = average[2]/mcs;
        AM2(counter) = average[3]/mcs;
        AabsM(counter) = average[4]/mcs;

        Evariance(counter) = (AE2(counter)- AE(counter)*AE(counter))/n_spins/n_spins;
        Mvariance(counter) = (AM2(counter) - AM(counter)*AM(counter))/n_spins/n_spins;
        M2variance(counter) = (AM2(counter) - AabsM(counter)*AabsM(counter))/n_spins/n_spins;

        HeatCv(counter) = Evariance(counter)/t(counter)/t(counter);
        Susceptibility(counter) = M2variance(counter)/t(counter);
//        Mcs(counter) = mcs;
        counter ++; //

    }
    AE.save("AE.txt", raw_ascii);

    cout << "Magnetization: " << average[4]/mcs << endl;
    cout << "HeatCv: " << Evariance/temperature/temperature << endl;
    cout << "Susceptibility: " << M2variance/temperature << endl;

    cout << "Lattice size: " << n_spins << endl;
    cout << "Mcs: " << mcs << endl;
    cout << "Average: " << *average/mcs << endl;
    cout << "Temp " << temperature << endl;  // temperature after loop: temperature+=temp_step

    return 0;

}


void output(int n_spins, int mcs, double temperature, double *average)


{

    // .txt file start
    ofstream writer( "Prosjekt 4d.txt" , ios::app ) ;
    //   Read to file
    writer << "RESULTS Prosjekt 4b:" << endl;
    writer << "Lattice size: " << n_spins << endl;
    writer << "Mcs: " << mcs << endl;
        writer << "Temperature: " << temperature << endl;
        writer << "Average: " << *average << endl;
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
void initialize(long& idum, int n_spins, double temperature, imat &spin_matrix,
                double& E, double& M)
{

    //     setup spin matrix and intial magnetization
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            //spin_matrix(y,x) = 1; // ORDER spin orientation for the ground state
            if (ran0(&idum)>0.5){
                spin_matrix(y,x) = 1;
            } else {
                spin_matrix(y,x) = -1;
            } // DISORDER spin orientation for the ground state
            M +=  (double) spin_matrix(y,x);
        }
    }

    // setup initial energy
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            E -=  (double) spin_matrix(y,x)*
                    (spin_matrix(periodic(y,n_spins,-1),x) +
                     spin_matrix(y,periodic(x,n_spins,-1)));
        }
    }
}// end function initialise


void Metropolis(int n_spins, long& idum, imat &spin_matrix, double& E, double&M, double *w)
{
    // loop over all spins
    // Step 1:Initial state: position random in lattice
    for(int i=0; i<n_spins*n_spins; i++) {
        int ix = (int) (ran0(&idum)*(double)n_spins);
        int iy = (int) (ran0(&idum)*(double)n_spins);

        int deltaE =  2*spin_matrix(iy,ix)*
                (spin_matrix(iy,periodic(ix,n_spins,-1))+
                 spin_matrix(periodic(iy,n_spins,-1),ix) +
                 spin_matrix(iy,periodic(ix,n_spins,1)) +
                 spin_matrix(periodic(iy,n_spins,1),ix));

        if ( ran0(&idum) <= w[deltaE+8] ) {
            spin_matrix(iy,ix) *= -1;  // flip one spin and accept new spin config
            M += (double) 2*spin_matrix(iy,ix);
            E += (double) deltaE;
        }
    }
} // end of Metropolis sampling over spins
