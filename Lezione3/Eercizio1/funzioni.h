#ifndef __funzioni__
#define __funzioni__

using namespace std;

// Funzione per calcolare l'errore statistico: prende in input due vettori di medie 
// e medie quadrate e il numero di blocchi, restituisce l'errore statistico
double error(vector<double>& av , vector<double>& av2, int n);

// Funzione per calcolare il chi-quadrato: Confronta i valori osservati 
// con quelli attesi e restituisce il valore del chi-quadrato
double compute_chi_squared(int expected, vector<int>& n_i);

// Funzione che verifica se una linea, lanciata a un angolo angolare, interseca una linea di distanza d
bool Intersect(double l, double d, double ang, Random& rnd);

// Funzione per calcolare la distanza quadratica
double distanza2(const vector<double>& pos);

//
void walker(vector<double>& posizione, double step, Random& rnd);

//
void walker_cont(vector<double>& posizione, double step, Random& rnd);

//
void initiating_rnd(Random& rnd);

#endif //__funzioni__