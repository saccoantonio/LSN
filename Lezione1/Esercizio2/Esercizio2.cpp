#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "random.h"
#include "funzioni.h"

using namespace std;

int main (){

    Random rnd;  initiating_rnd(rnd);
    
    unsigned int n = int(pow(10, 4)); // Numero di esperimenti
    vector<int> N_vector = {1, 2, 10, 100}; // Valori di N per il teorema del limite centrale
    double lambda = 1, mean = 0, gamma = 1; // Parametri per distribuzioni esponenziali e Lorentziane
    vector<double> st, exp, lor; // Vettori per salvare i risultati delle distribuzioni

    // Generazione delle distribuzioni per diversi valori di N
    for(unsigned int k = 0; k < N_vector.size(); k++){
        int N = N_vector[k];
        for(unsigned int i = 0; i < n; i++){
            double sum_s = 0, sum_e = 0, sum_l = 0;
            for(int j = 0; j < N; j++){
                sum_s += rnd.Rannyu(); // Distribuzione uniforme
                sum_e += rnd.Exponential(lambda); // Distribuzione esponenziale
                sum_l += rnd.Lorentz(mean, gamma); // Distribuzione Lorentziana
            }
            // Salvataggio delle medie per ogni distribuzione
            st.push_back(sum_s / N);
            exp.push_back(sum_e / N);
            lor.push_back(sum_l / N);
        }
    }
    
    // Scrittura dei dati generati su file
    ofstream distr;
    distr.open("../Files/distr.dat");
    if (!distr) {
        cerr << "Errore nell'apertura del file distr.dat" << endl;
        exit(1);
    }
    distr << "# unif(N=1)    exp(n=1)    lor(N=1)    unif(N=2)    exp(n=2)    lor(N=2)    unif(N=10)    exp(n=10)    lor(N=10)    unif(N=100)    exp(n=100)    lor(N=100)" << endl;
    for(unsigned int i = 0; i < n; i++){
        for(int j = 0; j < 4; j++){
            distr << st[i + j * n] << setw(15) << exp[i + j * n] << setw(15) << lor[i + j * n] << setw(15);
        }
        distr << endl;
    }
    distr.close();
    
    // Salvataggio del seed per riproducibilità
    rnd.SaveSeed();

    return 0;
}
