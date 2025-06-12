#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "random.h"
#include "funzioni.h"

using namespace std; 

double f(double x){
    return (M_PI/2)*cos(x*M_PI/2);
}
double p(double x){
    return 2*(1-x);
}

int main (int argc, char* argv[]){

    if (argc != 2) {
        cerr << "Errore: fornire un numero (1) o (2):" << '\n' 
        << "1) Metodo della Media" << '\n' 
        << "2) Importance Sampling" << endl;
        return 1;
    }
    int modalità = atoi(argv[1]);
    if (modalità == 1) {
        cerr << "Hai selezionato: 1) Metodo della Media" << endl;
    } else if(modalità == 2) {
        cerr << "Hai selezionato: 2) Importance Sampling" << endl;
    } else {
        cerr << "Errore: fornire un numero (1) o (2):" << '\n' 
        << "1) Metodo della Media" << '\n' 
        << "2) Importance Sampling" << endl;
        return 1;
    }
    
    Random rnd;  initiating_rnd(rnd);    

    int M = pow(10,4), L=100, integrazioni=100;
    int N = M/L; //numero di blocchi
    vector<double> sum, sum2;
    vector<double> var, var2; 

    if(modalità==1){
        for(int i=0; i<N; i++){
            double suma=0;
            double somma=0;
            for(int j=0; j<L; j++){
                for(int k=0; k<integrazioni; k++){
                    double t=rnd.Rannyu();
                    suma+=f(t)/integrazioni;
                    somma+=pow(f(t) ,2)/integrazioni;
                }
            }
            sum.push_back(suma/L);
            sum2.push_back(pow(suma/L,2));
            var.push_back(somma/L - pow(suma/L, 2));
            var2.push_back(pow(var[i], 2));
        }
    }else{
        for(int i=0; i<N; i++){
            double suma=0;
            double somma=0;
            for(int j=0; j<L; j++){
                for(int k=0; k<integrazioni; k++){
                    double t=rnd.cumulativa();
                    suma+=f(t)/p(t)/integrazioni;
                    somma+=pow(f(t)/p(t) ,2)/integrazioni;
                }
            }
            sum.push_back(suma/L);
            sum2.push_back(pow(suma/L,2));
            var.push_back(somma/L - pow(suma/L, 2));
            var2.push_back(pow(var[i], 2));
        }
    }

    // Calcolo delle medie progressive e degli errori
    vector<double> var_prog(N, 0.0), var2_prog(N, 0.0), std_prog(N, 0.0);
    vector<double> sum_prog(N, 0.0), sum2_prog(N, 0.0), err_prog(N, 0.0);  
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i + 1; j++){
            sum_prog[i] += sum[j];
            sum2_prog[i] += sum2[j];
            var_prog[i] += var[j];
            var2_prog[i] += var2[j];
        }
        sum_prog[i] /= (i + 1);
        sum2_prog[i] /= (i + 1);
        err_prog[i] = error(sum_prog, sum2_prog, i);
        var_prog[i] /= (i + 1);
        var2_prog[i] /= (i + 1);
        std_prog[i] = error(var_prog, var2_prog, i);
    }
    

    // Salvataggio dei risultati su file
    ofstream out;
    if(modalità==1){
        out.open("../Files/media.dat");
        cout << "file output -> media.dat: ";
    }else{
        out.open("../Files/sampling.dat");
        cout << "file output -> sampling.dat: ";
    }
    out << "# av /tab/ err /tab/ var /tab/ std" <<endl;
    for(int i = 0; i < N; i++){
        out << sum_prog[i] << "\t" << err_prog[i] << "\t" 
            << var_prog[i] << "\t" << std_prog[i] << endl;
    }
    out.close();

    rnd.SaveSeed();
    return 0;
}