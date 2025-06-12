#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "random.h"
#include "funzioni.h"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Errore: fornire un numero (1) o (2):" << '\n' 
             << "1) Cubic lattice random walk" << '\n' 
             << "2) Continuum random walk" << endl;
        return 1;
    }
    
    int modalità = atoi(argv[1]);
    if (modalità == 1) {
        cerr << "Hai selezionato: 1) Cubic lattice random walk" << endl;
    } else if (modalità == 2) {
        cerr << "Hai selezionato: 2) Continuum random walk" << endl;
    } else {
        cerr << "Errore: scegliere (1) o (2)" << endl;
        return 1;
    }

    Random rnd;  initiating_rnd(rnd);

    double step = 1.0;
    int passi = 100, lanci = 100000, blocchi = 100;
    int iter = lanci / blocchi;
    
    vector<double> sum(passi, 0.0), sum2(passi, 0.0), err(passi, 0.0);

    if (modalità == 1) {
        for (int k = 0; k < blocchi; k++) {
            vector<double> r2(passi, 0.0);
            for (int j = 0; j < iter; j++) {
                vector<double> posizione(3, 0.0);
                for (int i = 0; i < passi; i++) {
                    walker(posizione, step, rnd);
                    r2[i] += distanza2(posizione);  // Accumula |r_N|²
                }
            }
            for (int i = 0; i < passi; i++) {
                r2[i] /= iter;       // Normalizza su tutte le iterazioni del blocco
                sum[i] += r2[i];     // Accumula media nei blocchi
                sum2[i] += r2[i] * r2[i]; // Accumula il quadrato della media
            }
        }
    }else{
        for (int k = 0; k < blocchi; k++) {
            vector<double> r2(passi, 0.0);
            for (int j = 0; j < iter; j++) {
                vector<double> posizione(3, 0.0);
                for (int i = 0; i < passi; i++) {
                    walker_cont(posizione, step, rnd);
                    r2[i] += distanza2(posizione);  // Accumula |r_N|²
                }
            }
            for (int i = 0; i < passi; i++) {
                r2[i] /= iter;       // Normalizza su tutte le iterazioni del blocco
                sum[i] += r2[i];     // Accumula media nei blocchi
                sum2[i] += r2[i] * r2[i]; // Accumula il quadrato della media
            }
        }
    }

    ofstream out;
    if (modalità == 1) {
        out.open("../Files/cubic_lattice.dat");
        cout << "file output -> cubic_lattice.dat" << endl;
    } else {
        out.open("../Files/continuum.dat");
        cout << "file output -> continuum.dat" << endl;
    }
    out << "# √⟨|r_N|²⟩    err" << endl;
    for (int i = 0; i < passi; i++) {
        sum[i] /= blocchi;
        sum[i]=sqrt(sum[i]);
        sum2[i] /= blocchi;
        sum2[i]=sqrt(sum2[i]);
        err[i]=sqrt((sum2[i] - pow(sum[i], 2)) / 100);

        out << sum[i] << "\t" << err[i] << endl;
    }

    out.close();
    rnd.SaveSeed();
    return 0;
}
