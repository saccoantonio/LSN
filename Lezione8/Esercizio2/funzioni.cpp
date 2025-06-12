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

// Funzione per calcolare l'errore statistico: prende in input due vettori di medie 
// e medie quadrate e il numero di blocchi, restituisce l'errore statistico
double error(vector<double>& av, vector<double>& av2, int n){
    if(n==0){
        return 0;
    } else {
        return sqrt((av2[n] - pow(av[n], 2)) / n);
    }
}
double error_double( double AV, double AV2 , int n ){
    // funzione che mi restituisce la deviazione standard di n elementi
    if(n == 0){
        return 0;
    }
    return sqrt( (AV2 - AV*AV)/n );
    // oss: sarebbe diviso il numero di elementi - 1, ma io divido per n poiché l'indice dell'array parte da 0
}

// Funzione per calcolare il chi-quadrato: Confronta i valori osservati 
// con quelli attesi e restituisce il valore del chi-quadrato
double compute_chi_squared(int expected, vector<int>& n_i){
    double chi2 = 0.0;
    for (unsigned int i = 0; i < n_i.size(); i++) {
        chi2 += pow(n_i[i] - expected, 2) / expected;
    }
    return chi2;
}

// Funzione che verifica se una linea, lanciata a un angolo angolare, interseca una linea di distanza d
bool Intersect(double l, double d, double ang, Random& rnd){
    double y = (l/2) * sin(ang); // Calcola la coordinata y del punto di intersezione
    double x = rnd.Rannyu() * d; // Genera un punto casuale lungo la distanza d
    return y >= min(x, d - x); // Se la coordinata y è >= alla distanza tra x e d ==> intersezione avviene
}

// Funzione per calcolare la distanza quadratica
double distanza2(const vector<double>& pos) {
    double r2 = 0;
    for (unsigned int i = 0; i < pos.size(); i++) {
        r2 += pow(pos[i], 2);
    }
    return r2;
}

void walker(vector<double>& posizione, double step, Random& rnd){
    double num = rnd.Rannyu() * 3;
    if (num < 1) {
        if (rnd.Rannyu() < 0.5) {
            posizione[0] += step;
        } else {
            posizione[0] -= step;
        }
    } else if (num < 2) {
        if (rnd.Rannyu() < 0.5) {
            posizione[1] += step;
        } else {
            posizione[1] -= step;
        }
    } else {
        if (rnd.Rannyu() < 0.5) {
            posizione[2] += step;
        } else {
            posizione[2] -= step;
        }
    }
    return;
}

void walker_cont(vector<double>& posizione, double step, Random& rnd){
    double phi = rnd.Rannyu()*2*M_PI;
    double theta = rnd.Rannyu()*M_PI;
    posizione[0]+= step*sin(theta)*cos(phi);
    posizione[1]+= step*sin(theta)*sin(phi);
    posizione[2]+= step*cos(theta);    
    return;
}

void initiating_rnd(Random& rnd){
    int seed[4];
    int p1, p2;
    
    // Lettura dei numeri primi per il generatore di numeri casuali
    ifstream Primes("Primes");
    if (Primes.is_open()){
       Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
 
    // Lettura del seed iniziale
    ifstream input("seed.in");
    string property;
    if (input.is_open()){
       while (!input.eof()){
          input >> property;
          if (property == "RANDOMSEED"){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             rnd.SetRandom(seed, p1, p2);
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

double psi(double x, double mu, double sigma){
    return exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + exp(-(x+mu)*(x+mu)/(2*sigma*sigma));
}

double Hpsi(double x , double mu, double sigma){
    if (psi(x,mu,sigma) < 1e-10) return 1e6; // evita NaN
    return -0.5 *(  -1/(sigma*sigma) *  exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) -1/(sigma*sigma) *  exp(-(x+mu)*(x+mu)/(2*sigma*sigma)) + 1/(pow(sigma,4)) * (x-mu)*(x-mu) * exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + 1/(pow(sigma,4)) * (x+mu)*(x+mu) * exp(-(x+mu)*(x+mu)/(2*sigma*sigma)) )/psi(x, mu, sigma) + pow(x,4) - 5./2. * pow(x,2);
}
bool H_metro(double x , double y, double mu, double sigma, Random&  rnd){ // Metropolis algorithm
    bool decision = false;
    double psi_x = psi(x , mu , sigma);
    double psi_y = psi(y , mu , sigma);
    double acceptance = min( 1. , pow(psi_y/psi_x ,2) );
    if(rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
    return decision;
}

vector<double> H_measurement(double x, double mu, double sigma, Random&  rnd, int n_blocks, int n_steps, double& delta){

    int tune_steps = 1000; // numero di passi per il tuning del delta
    find_delta(x, mu, sigma, delta, tune_steps, rnd);

    double ave_blocco=0.0, ave2_blocco=0.0, prog_sum=0.0, prog_sum2=0.0;
    for(int i=0; i <n_blocks; i++){ 
        double stima = 0;

        for(int j=0; j < n_steps; j++){ 
            double y = x + rnd.Rannyu(-1.0,1.0) * delta;
            if(H_metro(x, y, mu, sigma, rnd)){
                x = y;
            }
            stima += Hpsi(x, mu, sigma);
        }
        prog_sum= prog_sum + stima/n_steps;
        prog_sum2 =prog_sum2 + stima/n_steps * stima/n_steps;
        ave_blocco = (prog_sum) / (i+1);
        ave2_blocco = ( prog_sum2 ) / (i+1);
    }
    vector <double> H;
    H.push_back(ave_blocco);
    H.push_back(error_double(ave_blocco, ave2_blocco, n_blocks-1));
    return H; 
}

void find_delta(double x0, double mu, double sigma, double& delta, int steps, Random&  rnd){
    
    double exp_acc=0.5; 
    double acc = 0;
    double x = x0;
    double tol = 0.01; 

    int max_iter = 1000;
    int iter = 0;
    while (iter < max_iter) {
        int accepted = 0;
        for (int i = 0; i < steps; ++i) {
            double y = x + rnd.Rannyu(-1., 1.) * delta;
            if (H_metro(x, y, mu, sigma, rnd)) {
                x = y;
                accepted++;
            }
        }
        acc = static_cast<double>(accepted) / steps;
        if (fabs(acc - exp_acc) < tol) {
            break;
        }
        if (acc < exp_acc)
            delta *= 0.9; // riduci se accettanza troppo bassa
        else
            delta *= 1.1; // aumenta se troppo alta

        iter++;
    }
}

bool params_metro(Random& rnd, double H_old , double H_new , double beta ){
    bool decision = false;
    double prob = exp(-beta * ( H_new - H_old ) );
    if(rnd.Rannyu() < prob ) decision = true;
    return decision;
}


void Execute_8_1(double mu, double sigma, double delta, int n_blocks, int n_steps){

    cout << "mu=" << mu << "  sigma=" << sigma << "  delta=" << delta << "  n_blocks=" << n_blocks << "  n_steps=" << n_steps << endl; 

    Random rand; initiating_rnd(rand);
    
    ofstream acc_file("../Files_8.2/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
    acc_file << "#   N_BLOCK:  ACCEPTANCE:" << endl;
    ofstream energy_file("../Files_8.2/energy.dat");
    energy_file << "#     BLOCK:  ACTUAL_E:     E_AVE:      ERROR:" << endl;
    ofstream psi_file("../Files_8.2/psi_sampled.dat");
    psi_file << "# Sampled positions:" << endl;

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
            double y = x +  rand.Rannyu(-1.0,1.0) * delta;
            if(H_metro(x,y,mu, sigma,rand)){
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
}

