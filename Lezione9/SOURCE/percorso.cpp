#include "percorso.h"

using namespace std;

vector<vector<double>> percorso::matrice_distanze;

percorso :: percorso(int n_citta, Random& r): n_city(n_citta), path(n_citta), rnd(&r){
    for(int i=0; i<n_city; i++){
        path[i]=i+1;
    }
}

percorso :: ~percorso(){}

bool percorso :: bounds_control(){ // verifica che path[0]=1 e che ogni città è percorsa una sola volta
    bool ripetizioni = false;
    vector<bool> visto(n_city, false);
    for (int x : path) {
        if (visto[x]) ripetizioni = true; // già visto
        visto[x] = true;
    }
    bool control;
    if(path[0]==1 && ripetizioni == false){
        control = true ;
    } else {control = false;}
    return control;
}

void percorso :: set_city(int position, int city){
    path[position] = city;
};

int percorso :: get_city(int position){
    return path[position];
};

int percorso :: get_n_city(){
    return n_city;
}

int percorso :: pbc(int posto) {
    return (posto % n_city + n_city) % n_city;
};

void percorso :: pair_swap(){
    int posto1 = floor(rnd->Rannyu()*(n_city-1)) + 1;
    int posto2 = floor(rnd->Rannyu()*(n_city-1)) + 1;
    swap(path[posto1], path[posto2]);
    if(bounds_control() == false) cout << "vincoli non rispettati" << endl;
};

void percorso :: shift(){
    vector<int> temp(n_city - 1);
    for (int i = 1; i<n_city; i++) {
        temp[i - 1] = path[i];
    }
    int m = floor(rnd->Rannyu() * (n_city - 1));

    int mod = n_city-1;
    for (int i = 1; i<n_city; i++) {
        int shifted_index = ((i - 1 - m) % mod + mod) % mod;;
        path[i] = temp[shifted_index];
    }
    if(bounds_control() == false) cout << "vincoli non rispettati" << endl;
};

void percorso :: block_inversion() { //all'inizio nn avevo messo le pbc
    int start = floor(rnd->Rannyu() * (n_city - 1)) + 1; // da 1 a n_city-1
    int end = floor(rnd->Rannyu() * (n_city - 1)) + 1;   // da 1 a n_city-1
    if (start == end) end = pbc(end + 4); 
    if (end == 0) end = 1;

    int len = (end - start + (n_city - 1)) % (n_city - 1) + 1;
    for (int i = 0; i < len / 2; ++i) { // se il blocco è dispari fa niente
        int id1 = pbc(start + i);
        int id2 = pbc(end - i);
        if (id1 == 0) id1 = 1;
        if (id2 == 0) id2 = 1;
        swap(path[id1], path[id2]);
    }
    if (bounds_control() == false) cout << "vincoli non rispettati" << endl;
}


void percorso :: block_permutation() {
    int start = floor(rnd->Rannyu() * (n_city - 1)) + 1;
    int end   = floor(rnd->Rannyu() * (n_city - 1)) + 1;
    if (start == end) end = pbc(end + 4); 
    if (end == 0) end = 1;

    int m = floor(rnd->Rannyu() * n_city * 2);
    for (int i=0; i<m; i++) {
    int i1 = pbc(floor(rnd->Rannyu() * abs(end - start + 1)) + start);
    int i2 = pbc(floor(rnd->Rannyu() * abs(end - start + 1)) + start);
    if (i1==0) i1++;
    if (i2==0) i2++;
    swap(path[i1], path[i2]);
    }
    if(bounds_control() == false) cout << "vincoli non rispettati" << endl;
};

void percorso :: multiple_mutation(){ // questa non funziona troppo bene
    if(rnd->Rannyu() <= 0.1) pair_swap();
    if(rnd->Rannyu() <= 0.05) shift();
    if(rnd->Rannyu() <= 0.05) block_inversion();
    if(rnd->Rannyu() <= 0.05) block_permutation();    
};

void percorso :: single_mutation(){
    if(rnd->Rannyu() <= 0.2){
        double mutaz = rnd->Rannyu();
        if (mutaz<0.25) pair_swap();
        else if (mutaz<0.5) shift();
        else if (mutaz<0.75) block_inversion();
        else block_permutation();
    }
}

double percorso::get_lenght() {
    double lunghezza = 0.0;
    for (int i = 0; i < n_city; i++) {
        int index1 = path[i] - 1;
        int index2 = path[pbc(i + 1)] - 1;
        lunghezza += (matrice_distanze)[index1][index2];
    }
    return lunghezza;
}

Random* percorso::get_rnd() {
    return rnd;
}

vector<int> percorso :: get_path(){
    return path;
};