#include "generazione.h"


generazione :: generazione(int n_citta, int n_ind, int r): n_city(n_citta), n_individui(n_ind), rank(r){ 
    initiating_rnd_wline(rnd, rank); 
    for (int i=0; i<n_individui; i++) {
        percorso p(n_city, rnd);
        individui.push_back(p);
    }
}

generazione :: ~generazione(){};

int generazione :: get_n_individui(){
    return n_individui;
};

int generazione :: get_n_city(){
    return individui[0].get_n_city();
}

vector<percorso> generazione :: get_all_paths(){
    return individui;
};

percorso generazione :: get_path(int posto){
    return individui[posto];
}

void generazione :: set_path(int posizione, percorso path){
    individui[posizione] = path;
};

percorso generazione :: get_best(){
    ordina();
    return individui[0];
};

void generazione :: ordina() {
    sort(individui.begin(), individui.end(), 
        // a quanto pare non posso usare una funz definita da me ma solo una lambda
        [=](percorso& a, percorso& b) {
            return a.get_lenght() < b.get_lenght();
        }
    );
}

void generazione::order_crossover() {

    ordina();

    vector<percorso> figli(n_individui, percorso(n_city, rnd));

    int figlio_index = 0;
    int accoppiamenti = floor(n_individui / 2);

    for (int k = 0; k < accoppiamenti; k++) {
    int index_mother = floor(n_individui * pow(rnd.Rannyu(), 2));
    int index_father = floor(n_individui * pow(rnd.Rannyu(), 2));

    percorso mother = get_path(index_mother);
    percorso father = get_path(index_father);

    if (rnd.Rannyu() >= 0.9) { // lascio la scelta ai genitori di accoppiarsi o meno
        figli[figlio_index++] = mother;
        figli[figlio_index++] = father;
        continue;
    }

    int n_city = mother.get_n_city();
    int cut = floor(rnd.Rannyu(1, n_city - 2));  // cut tra 1 e n_city-2

    percorso child1(n_city, rnd);
    percorso child2(n_city, rnd);

    for (int i = 0; i < cut; i++) {
        child1[i] = mother[i];
        child2[i] = father[i];
    }

    vector<bool> citta_viste1(n_city, false);
    vector<bool> citta_viste2(n_city, false);
    for (int i = 0; i < cut; i++) {
        citta_viste1[child1[i] - 1] = true;
        citta_viste2[child2[i] - 1] = true;
    }

    vector<int> temp1, temp2;
    for (int i = 0; i < n_city; i++) {
        if (!citta_viste1[father[i] - 1]) temp1.push_back(father[i]);
        if (!citta_viste2[mother[i] - 1]) temp2.push_back(mother[i]);
    }

    for (int i = 0; i < n_city - cut; i++) {
        child1[cut + i] = temp1[i];
        child2[cut + i] = temp2[i];
    }

    figli[figlio_index++] = child1;
    figli[figlio_index++] = child2;
}

    // Se n_individui è dispari, copia uno dei migliori per riempire lo slot mancante
    if (n_individui % 2 == 1) {
        figli[n_individui - 1] = get_path(0);  // copia il migliore
    }

    // Controllo di vincoli opzionale
    for (int i = 0; i < n_individui; i++) {
        if (!figli[i].bounds_control()) {
            cout << "Vincoli non rispettati nel figlio " << i << endl;
        }
    }

    individui = figli;
}

void generazione :: global_mutation(){
    for(int i=0; i<n_individui; i++){
        individui[i].single_mutation();
    }
};

void generazione::many_crossovers(int accoppiamenti) { 
    for(int i=0; i<accoppiamenti; i++){
        ordina();
        order_crossover();
    }    
}

double generazione :: average_distance(){
    double somma=0;
    for (int i=0; i<floor(n_individui/2); i++){
        somma+=individui[i].get_lenght();
    }
    return somma/floor(n_individui/2);
}