#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "random.h"
#include "funzioni.h"

using namespace std; 

int main (){

    Random rnd;  initiating_rnd(rnd);

    // ~~~~~~~~~~~ Punto 1 + 2 ~~~~~~~~~~~ //
    
    int N = 100, L = 1000; // Numero totale di dati, blocchi e grandezza blocco
    vector<double> av, av2;
    vector<double> sum_prog(N, 0.0), sum2_prog(N, 0.0), err_prog(N, 0.0);
    vector<double> var, var2; 
    vector<double> var_prog(N, 0.0), var2_prog(N, 0.0), std_prog(N, 0.0);

    // Ciclo per la media e la varianza
    for(int i = 0; i < N; i++){
        double sum = 0, somma = 0;
        for(int j = 0; j < L; j++){
            sum += rnd.Rannyu(); // Genera numero casuale uniforme
            somma += pow(rnd.Rannyu() - 0.5, 2); // Varianza
        }
        av.push_back(sum / L);
        av2.push_back(pow(av[i], 2));
        var.push_back(somma / L);
        var2.push_back(pow(var[i], 2));
    }
    
    // Calcolo delle medie progressive e degli errori
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i + 1; j++){
            sum_prog[i] += av[j];
            sum2_prog[i] += av2[j];
            var_prog[i] += var[j];
            var2_prog[i] += var2[j];
        }
        sum_prog[i] /= (i + 1);
        sum2_prog[i] /= (i + 1);
        err_prog[i] = error(sum_prog, sum2_prog, i);
        var_prog[i] /= (i + 1);
        var2_prog[i] /= (i + 1);
        std_prog[i] = error(var_prog, var2_prog, i);
    }

    // Salvataggio dei risultati su file
    ofstream out;
    out.open("../Files/data.dat");
    out << "# av /tab/ err /tab/ var /tab/ std" <<endl;
    for(int i = 0; i < N; i++){
        out << sum_prog[i] << "\t" << err_prog[i] << "\t" 
            << var_prog[i] << "\t" << std_prog[i] << "\t" 
            << endl;
    }
    out.close();
    
    // ~~~~~~~~~~~~~~ Punto 3 ~~~~~~~~~~~~~~ //
    
    int M = 100; // Numero di intervalli per il chi-quadrato
    int n = 10000, exp = n / M; // Numero di estrazioni e valore atteso per intervallo
    int tests = 1000;
    ofstream chi;
    chi.open("../Files/chi2.dat");
   
    // Calcolo del chi-quadrato
    for(int j = 0; j < tests; j++){
        vector<int> n_i(M, 0); // Conta il numero di eventi in ogni intervallo
        for(int i = 0; i < n; i++){
            int interval = floor(rnd.Rannyu() * M);
            n_i[interval]++;
        }
        chi << compute_chi_squared(exp, n_i) << endl;
    }
    chi.close();
    
    rnd.SaveSeed(); // Salvataggio del seed per riproducibilità
    return 0;
}
