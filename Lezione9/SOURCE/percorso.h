#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include "random.h"
#include "funzioni.h"
#ifndef __percorso__
#define __percorso__

using namespace std;

class percorso {
private:
    int n_city;
    vector<int> path;
    Random* rnd;
    static vector<vector<double>> matrice_distanze; // grazie chat <3
public:
    percorso(int n_citta, Random& r);
    ~percorso();
    static void set_distance_matrix(const vector<vector<double>>& d){
        matrice_distanze = d;
    };
    bool bounds_control();
    void set_city(int position, int city);
    int get_city(int position);
    int get_n_city();
    int pbc(int posto);  
    void pair_swap();
    void shift();
    void block_inversion();
    void block_permutation();
    void multiple_mutation();
    void single_mutation();
    double get_lenght();
    int& operator[](int i) { return path[i]; } //inserito per poter usare la scrittura percorso[i]
    const int& operator[](int i) const { return path[i]; } 
    percorso& operator=(const percorso&) = default;
    Random* get_rnd();
    vector<int> get_path();
};

#endif // __percorso__