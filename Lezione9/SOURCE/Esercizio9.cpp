#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include "random.h"
#include "percorso.h"
#include "generazione.h"
#include "funzioni.h"

using namespace std;

int main() {
    ifstream input("../sim_parameters.dat");
    if (!input.is_open()) {
        cerr << "Errore: impossibile aprire il file." << endl;
        return 1;
    }

    string key;
    int n_individui = 0, n_citta = 0, n_generazioni = 0;
    string nomefilein;
    ofstream L1, best_path;

    while (input >> key) {
        if (key == "NCITTA") {
            input >> n_citta;
        } else if (key == "NINDIVIDUI") {
            input >> n_individui;
        } else if (key == "NGENERATIONS") {
            input >> n_generazioni;
        } else if (key == "NOMEFILEIN") {
            input >> nomefilein;
            nomefilein = "../CONFIGS/" + nomefilein;

            if (nomefilein == "../CONFIGS/circular_path.txt") {
                L1.open("../OUTPUT/L1_circ.dat");
                best_path.open("../OUTPUT/best_circ.dat");
            } else if (nomefilein == "../CONFIGS/random_path.txt") {
                L1.open("../OUTPUT/L1_rand.dat");
                best_path.open("../OUTPUT/best_rand.dat");
            } else if (nomefilein == "../CONFIGS/cap_prov_ita.dat") {
                L1.open("../OUTPUT/L1_cap.dat");
                best_path.open("../OUTPUT/best_cap.dat");
            }else {
                cout << "Errore: nome file non riconosciuto." << endl;
                exit(1);
            }

            L1 << "# GEN: \t LUNGHEZZA DEL MIGLIORE: \t LUNGHEZZA MEDIA:" << endl;
            best_path << "# Percorso migliore trovato" << endl;
        }
    }
    input.close();

    ifstream posizioni(nomefilein);
    if (!posizioni.is_open()) {
        cout << "Errore: impossibile aprire il file." << endl;
        return -1;
    }

    generazione pop(n_citta, n_individui, 0);
    vector<vector<double>> citta(n_citta, vector<double>(2, 0.0));

    if(nomefilein == "../CONFIGS/cap_prov_ita.dat"){ // inserito per confrontarmi con il programma in parallelo
        string line;
        for(int i=0; i<n_citta; i++){
            getline(posizioni, line);
            istringstream iss(line);
            double x, y;
            iss >> x >> y;
            citta[i][0] = x;
            citta[i][1] = y;
        }
    } else {
        string line;
        getline(posizioni, line); // Salta intestazione
        while (getline(posizioni, line)) {
            istringstream iss(line);
            int id; double x, y;
            if (iss >> id >> x >> y) {
                if (id >= 1 && id <= n_citta) {
                    citta[id - 1][0] = x;
                    citta[id - 1][1] = y;
                }
            }
        }
    }

    posizioni.close();
    vector<vector<double>> distance_matrix(n_citta, vector<double>(n_citta, 0.0));
    for (int i = 0; i < n_citta; ++i) {
        for (int j = 0; j < n_citta; ++j) {
            distance_matrix[i][j] = sqrt( pow(citta[i][0] - citta[j][0], 2) + pow(citta[i][1] - citta[j][1], 2) );
        }
    }
    percorso::set_distance_matrix(distance_matrix);

    for (int gen = 0; gen < n_generazioni; gen++) {
        pop.global_mutation();
        pop.order_crossover();
        L1 << gen + 1 << "\t" << pop.get_best().get_lenght() << "\t" << pop.average_distance() << endl;
    }

    percorso best = pop.get_best();
    for (int i = 0; i < n_citta; i++) {
        best_path << best[i] << endl;
    }

    best_path.close();
    L1.close();

    return 0;
}
