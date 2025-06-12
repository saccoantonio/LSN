#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <armadillo>
#include <algorithm>
#include "random.h"
#include "funzioni.h"

using namespace std;

int main (int argc, char *argv[]){

    // random number generator
    Random rnd; initiating_rnd(rnd);
    
    ofstream acc_file("../Files/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
    acc_file << "#   N_BLOCK:  ACCEPTANCE:" << endl;
    ofstream energy_file("../Files/energy.dat");
    energy_file << "#     BLOCK:  ACTUAL_E:     E_AVE:      ERROR:" << endl;
    ofstream psi_file("../Files/psi_sampled.dat");
    psi_file << "# Sampled positions:" << endl;

    // Simulation parameters
    ifstream input("sim_parameters.dat");
    if (!input.is_open()) {
        cerr << "Errore: impossibile aprire il file." << endl;
        return 1;
    }
    string key;
    int n_steps = 0, n_blocks;
    double delta = 0.0, mu = 0.0, sigma = 0.0;
    while (input >> key) {
        if (key == "NSTEPS") {
            input >> n_steps;
        } else if (key == "NBLOCK") {
            input >> n_blocks;
        } else if (key == "MU") {
            input >> mu;
        } else if (key == "SIGMA") {
            input >> sigma;
        }else if (key == "DELTA") {
            input >> delta;
        }
    }
    input.close();


    double ave_blocco=0.0; // quando inserito nel ciclo for, contiene la media delle medie dei primi (i+1) blocchi
    double ave2_blocco=0.0; // quando inserito nel ciclo for, la media dei quadrati delle medie dei primi (i+1) blocchi
    double prog_sum=0.0; // accumula la somma delle medie dei primi (i+1) blocchi
    double prog_sum2=0.0; // accumula la somma delle medie dei quadrati 

    double n_attemps=0.0;
    double n_accepted=0.0;
    double x = 0.0; // configurazione (posizione) iniziale

    vector<double> acc;

    for(int i=0; i <n_blocks; i++){ //loop over blocks
        double stima = 0;
        n_attemps = 0;
        n_accepted = 0;
        for(int j=0; j < n_steps; j++){ //loop over steps in a block
            // step Metropolis 
            n_attemps++;
            double y = x +  rnd.Rannyu(-1.0,1.0) * delta;
            if(metro(x,y,mu, sigma,rnd)){
                x = y;
                n_accepted++;
                psi_file << x << endl;
            }
            stima +=  Hpsi(x, mu, sigma);
        }
        acc_file << (i+1) << "\t" << double(n_accepted)/double(n_attemps) << endl;
        prog_sum= prog_sum + stima/n_steps;
        prog_sum2 =prog_sum2 + stima/n_steps * stima/n_steps;
        ave_blocco = (prog_sum) / (i+1);
        ave2_blocco = ( prog_sum2 ) / (i+1);
        energy_file << (i+1) << "\t" << stima/n_steps  << "\t" << (ave_blocco ) << "\t" << error_double(ave_blocco, ave2_blocco, i) << endl;

        acc.push_back(double(n_accepted)/double(n_attemps));
    }
    double somma = accumulate(acc.begin(), acc.end(), 0.0);
    cout << "Accettanza media = " << somma / acc.size() << " con delta = " << delta << endl;
    energy_file.close();
    acc_file.close();
    psi_file.close();

    return 0;
}
