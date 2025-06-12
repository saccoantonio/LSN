#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"
#include "funzioni.h"

using namespace std;

int main(){
    Random rnd;  initiating_rnd(rnd);
    
    // ~~~~~~~~~~~~~~~~~~~~ Calcolo di pi greco con il metodo di Buffon ~~~~~~~~~~~~~~~~~~~~ // 

    int M = 1000000, N = 100, L = M / N, l = 1;  // M = numero di lanci, N = numero di esperimenti
    double d = 3;  // Lunghezza della distanza tra le linee

    // Vettori per memorizzare i risultati
    vector<double> sum_pi(N), sum_pi2(N), err_pi(N);
    
    // Loop su ciascun esperimento
    for(int i = 0; i < N; i++){
        int Nhit = 0;  // Numero di intersezioni per esperimento
        for(int j = 0; j < L; j++){
            double theta = Angle(rnd)/2;  // Angolo casuale tra 0 e pi
            if(Intersect(l, d, theta, rnd)) {
                Nhit++;  // Conta il numero di intersezioni
            }
        }
        sum_pi[i] = (2.0 * l * L) / (d * Nhit);  // Calcola la stima di pi per ogni esperimento
        sum_pi2[i] = pow(sum_pi[i], 2);  // Per il calcolo dell'errore
    }

    // Scrittura dei risultati su file
    ofstream out("../Files/pi.dat");
    out << "# stima_pi  err" << endl;
    vector<double> sum_prog(N, 0.0), sum_prog2(N, 0.0);
    
    // Calcolo delle medie e dell'errore per ogni esperimento
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            sum_prog[i] += sum_pi[j];  // Somma delle stime di pi
            sum_prog2[i] += sum_pi2[j];  // Somma dei quadrati delle stime di pi
        }
        sum_prog[i] /= (i + 1);  // Media delle stime fino all'indice i
        sum_prog2[i] /= (i + 1);  // Media dei quadrati delle stime fino all'indice i
        err_pi[i] = error(sum_prog, sum_prog2, i); // Calcola l'errore 
        out << sum_prog[i] << "\t" << err_pi[i] << endl;  // Scrive i risultati nel file
    }
    out.close();
    return 0;
}