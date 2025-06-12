#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"
#include "funzioni.h"

using namespace std;

// Funzione per calcolare l'errore statistico: prende in input due vettori di medie 
// e medie quadrate e il numero di blocchi, restituisce l'errore statistico
double error(vector<double> av , vector<double> av2, int n){
    if(n==0){
        return 0;
    } else {
        return sqrt((av2[n] - pow(av[n], 2)) / n);
    }
}

// Funzione per calcolare il chi-quadrato: Confronta i valori osservati 
// con quelli attesi e restituisce il valore del chi-quadrato
double compute_chi_squared(int expected, vector<int> n_i){
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

//Funzione che estrae angolo da 0 a 2PI
double Angle(Random& rnd){
    double x=rnd.Rannyu(-1,1);
    double y=rnd.Rannyu(-1,1);
    double theta;
    if(y>=0){
        theta=acos(x/sqrt(x*x + y*y));
    }else{
        theta=2*M_PI - acos(x/sqrt(x*x + y*y));
    }
    return theta;
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