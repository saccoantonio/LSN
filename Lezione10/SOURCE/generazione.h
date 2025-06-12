#include "percorso.h"

#ifndef __generazione__
#define __generazione__

using namespace std;

// This class contains functions for generating random numbers using the RANNYU algorithm
class generazione {
private:
    int n_city;
    int n_individui;
    vector<percorso> individui;
    Random rnd;    
    int rank;
public:
    generazione(int n_citta, int n_ind, int r);
    ~generazione();
    void ordina();
    int get_n_individui();
    int get_n_city();
    vector<percorso> get_all_paths();
    percorso get_path(int posto);
    void set_path(int posizione, percorso path);
    percorso get_best();
    void order_crossover();
    void many_crossovers(int accoppiamenti);
    void global_mutation(); 
    percorso& operator[](int i) { return individui[i]; } //inserito per poter usare la scrittura individui[i] per recuperare easy un percorso
    double average_distance();
};

#endif // __generzione__
