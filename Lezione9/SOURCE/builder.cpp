#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <armadillo>
#include <algorithm>
#include "random.h"
#include "funzioni.h"

using namespace std;

int main() {

    Random rnd; initiating_rnd(rnd);
    int n_city = 34;
    vector<int> citta (n_city);
    citta[0] = 1;
    for(int k=1; k<n_city; k++){
        citta[k] = k+1;
    }

    for(int k=0; k<n_city*10; k++){
        int posto1 = floor(rnd.Rannyu()*(n_city));
        int posto2 = floor(rnd.Rannyu()*(n_city));
        if (citta[posto1]!=1 && citta[posto2]!=1) swap(citta[posto1], citta[posto2]);
    }

    ofstream circ ("../CONFIGS/circular_path.txt");
    circ << "#CITTA' \t X \t Y" << endl;
    circ << citta[0] << "\t" << 0.0 << "\t" << 1 << endl;
    for(int k=1; k<n_city; k++){
        double theta = rnd.Rannyu(0, 360);
        circ << citta[k] << "\t" << cos(theta) << "\t" << sin(theta) << endl;
    }
    circ.close();

    ofstream random ("../CONFIGS/random_path.txt");
    random << "#CITTA' \t X \t Y" << endl;
    for(int k=0; k<n_city; k++){
        random << citta[k] << "\t" << rnd.Rannyu(-1,1) << "\t" << rnd.Rannyu(-1,1) << endl;
    }
    random.close(); 

    return 0;
}