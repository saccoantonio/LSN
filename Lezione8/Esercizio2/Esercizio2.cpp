#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "funzioni.h"
#include "random.h"

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd;  initiating_rnd(rnd);
    double x0 = 0.0;  // punto di partenza per Metropolis

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

    // Valore iniziale dell'energia
    vector<double> H = H_measurement(x0, mu, sigma, rnd, n_blocks, n_steps, delta);

    // File output
    ofstream par_file("../Files_8.2/evo_parameters.dat");
    par_file << "# T\t\tmu\t\tsigma\t\tH\t\terr" << endl;

    // Ciclo di simulated annealing
    for (double T = 2.0; T >= 0.01; T *= 0.99) {
        for (int n = 0; n < 100; n++) {
            double delta_mu = 0.5 * T;
            double delta_sigma = 0.5 * T;
            double beta = 1.0 / T;
            // Proposta per mu
            double mu_new = fabs(mu + rnd.Rannyu(-1, 1) * delta_mu);

            // Proposta per sigma con controllo
            double sigma_new;
            do {
                sigma_new = sigma + rnd.Rannyu(-1, 1) * delta_sigma;
            } while (sigma_new < 0.1);  // evitare sigma troppo piccole

            // Calcolo energia con nuovi parametri
            vector<double> H_new = H_measurement(x0, mu_new, sigma_new, rnd, n_blocks, n_steps, delta);

            // Accettazione Metropolis
            bool metro = params_metro(rnd, H[0], H_new[0], beta);
            if (metro) {
                mu = mu_new;
                sigma = sigma_new;
                H = H_new;
            }
            // Altrimenti tieni H corrente (già fatto)
        }
        // Scrivi parametri aggiornati e valore energia
        par_file << T << "\t\t" << mu << "\t\t" << sigma << "\t\t" << H[0] << "\t\t" << H[1] << endl;
        cout << T << endl;
    }
    par_file.close();

    Execute_8_1(mu, sigma, delta, n_blocks, n_steps);

    return 0;
}
