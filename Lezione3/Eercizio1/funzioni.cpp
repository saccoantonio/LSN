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
double error(vector<double>& av, vector<double>& av2, int n){
    if(n==0){
        return 0;
    } else {
        return sqrt((av2[n] - pow(av[n], 2)) / n);
    }
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