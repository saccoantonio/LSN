#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "random.h"
#include "funzioni.h"

using namespace std; 

int main (int argc, char* argv[]){

    if (argc != 2) {
        cerr << "Errore: fornire un numero (1) o (2):" << '\n' 
        << "1) Metodo diretto" << '\n' 
        << "2) Metodo discretizzato" << endl;
        return 1;
    }
    int modalità = atoi(argv[1]);
    if (modalità == 1) {
        cerr << "Hai selezionato: 1) Metodo diretto" << endl;
    } else if(modalità == 2) {
        cerr << "Hai selezionato: 2) Metodo discretizzato" << endl;
    } else {
        cerr << "Errore: fornire un numero (1) o (2):" << '\n' 
        << "1) Metodo diretto" << '\n' 
        << "2) Metodo discretizzato" << endl;
        return 1;
    }

    // ~~~~~~~~~~~ Econophysics ~~~~~~~~~~~ //
    
    Random rnd;  initiating_rnd(rnd);
    
    int blocchi=100, iter=1000, intervalli=100;
    double T=1.0, K=100, r=0.1, sigma=0.25, S0=100, t=T/intervalli;
    vector<double> av(blocchi, 0.0), av2(blocchi, 0.0);
    vector<double> med(blocchi, 0.0), med2(blocchi, 0.0);
    
    if (modalità == 1) { 
        for(int i = 0; i < blocchi; i++){
            double C=0, P=0;
            for(int j = 0; j < iter; j++){
                double Z = rnd.Gauss(0, 1);
                double S = S0*exp((r-pow(sigma,2)/2)*T+sigma*Z*sqrt(T));
                C += exp(-r*T)*max(0.0, S-K);
                P += exp(-r*T)*max(0.0, K-S);
            }
            av[i] = C/iter;
            av2[i] = pow(av[i], 2);
            med[i] = P/iter;
            med2[i] = pow(med[i], 2);
        }
    }else{
        for(int i = 0; i < blocchi; i++){
            double C=0, P=0, S=0;
            for(int j = 0; j < iter; j++){
                S=S0;
                for(int k=0; k<intervalli; k++){
                    double Z=rnd.Gauss(0, 1);
                    S *= exp((r-pow(sigma,2)/2)*t+sigma*Z*sqrt(t));
                }
                C += exp(-r*T)*max(0.0, S-K);
                P += exp(-r*T)*max(0.0, K-S);
            }
            av[i] = C/iter;
            av2[i] = pow(av[i], 2);
            med[i] = P/iter;
            med2[i] = pow(med[i], 2);
        }
    }

    // Calcolo delle medie progressive e degli errori
    vector<double> sum_prog(blocchi, 0.0), sum2_prog(blocchi, 0.0), err_prog(blocchi, 0.0);
    vector<double> tot_prog(blocchi, 0.0), tot2_prog(blocchi, 0.0), inc_prog(blocchi, 0.0);

    for(int i = 0; i < blocchi; i++){
        for(int j = 0; j < i + 1; j++){
            sum_prog[i] += av[j];
            sum2_prog[i] += av2[j];
            tot_prog[i] += med[j];
            tot2_prog[i] += med2[j];
        }
        sum_prog[i] /= (i + 1);
        sum2_prog[i] /= (i + 1);
        err_prog[i] = error(sum_prog, sum2_prog, i);

        tot_prog[i] /= (i + 1);
        tot2_prog[i] /= (i + 1);
        inc_prog[i] = error(tot_prog, tot2_prog, i);
    }
    
    // Salvataggio dei risultati su file
    ofstream out;
    if(modalità==1){
        out.open("../Files/diretto.dat");
        cout << "file output -> diretto.dat: ";
    }else{
        out.open("../Files/discretizzato.dat");
        cout << "file output -> discretizzato.dat: ";
    }
    out << "# C /tab/ C_err /tab/ P /tab/ P_err" << endl;
    for(int i = 0; i < blocchi; i++){
        out << sum_prog[i] << "\t" << err_prog[i] << "\t" << tot_prog[i] << "\t" << inc_prog[i] <<  endl;
    }
    out.close();
    
    rnd.SaveSeed(); // Salvataggio del seed per riproducibilità
    return 0;
}
