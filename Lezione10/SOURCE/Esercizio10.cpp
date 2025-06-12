#include "generazione.h"
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ifstream input("../sim_parameters.dat");
    if (!input.is_open()) {
        cerr << "Errore: impossibile aprire il file." << endl;
        MPI_Finalize();
        return 1;
    }

    string key;
    int n_individui = 0, n_citta = 0, n_generazioni = 0, n_migranti = 0;
    string nomefilein;
    ofstream L1, best_path;

    while (input >> key) {
        if (key == "NCITTA") {
            input >> n_citta;
        } else if (key == "NINDIVIDUI") {
            input >> n_individui;
        } else if (key == "NGENERATIONS") {
            input >> n_generazioni;
        } else if (key == "NMIGRANTI") {
            input >> n_migranti;
        } else if (key == "NOMEFILEIN") {
            input >> nomefilein;
            nomefilein = "../CONFIGS/" + nomefilein;
        }
    }
    input.close();

	for(int i=0; i<size; i++){
		if (rank == i) {
			L1.open("../OUTPUT/OUT"+to_string(rank)+"/L1.dat");
			best_path.open("../OUTPUT/OUT"+to_string(rank)+"/best.dat");
			L1 << "# GEN: \t LUNGHEZZA DEL MIGLIORE: \t LUNGHEZZA MEDIA:" << endl;
			best_path << "# Percorso migliore trovato dal core " + to_string(rank) << endl;
		}
	}

    ifstream posizioni(nomefilein);
    if (!posizioni.is_open()) {
        cout << "Errore: impossibile aprire il file." << endl;
        MPI_Finalize();
        return -1;
    }

    vector<vector<double>> citta(n_citta, vector<double>(2, 0.0));
    string line;
    for(int i=0; i<n_citta; i++){
        getline(posizioni, line);
        istringstream iss(line);
        double x, y;
        iss >> x >> y;
        citta[i][0] = x;
        citta[i][1] = y;
    }
    posizioni.close();

    vector<vector<double>> distance_matrix(n_citta, vector<double>(n_citta, 0.0));
    for (int i = 0; i < n_citta; ++i) {
        for (int j = 0; j < n_citta; ++j) {
            distance_matrix[i][j] = sqrt( pow(citta[i][0] - citta[j][0], 2) + pow(citta[i][1] - citta[j][1], 2) );
        }
    }

    percorso::set_distance_matrix(distance_matrix);
    
// ~=~-~o@o~-~=~<§][§>~=~-~o@o~-~=~<§][§>~=~-~o@o~-~=~<§][§>~=~-~o@o~-~=~<§][§>~=~-~o@o~-~=~<§][§>~=~-~o@o~-~=~<§][§>~=~-~o@o~-~=~ //

    generazione pop(n_citta, n_individui, rank);
    Random rand; initiating_rnd_wline(rand, rank);
    int rate_migrazione = 100;

    for (int gen = 0; gen < n_generazioni; gen++) {
        pop.global_mutation();
        pop.order_crossover();
        L1 << gen + 1 << "\t" << pop.get_best().get_lenght() << "\t" << pop.average_distance() << endl;

        // migrazione
        if (gen % rate_migrazione == 0 && size > 1) {
            for (int j = 0; j < n_migranti; j++) {

                int index = floor(n_individui * pow(rand.Rannyu(), 2));
				percorso path = pop.get_path(index);
				vector<int> buffer_out(n_citta);
				for (int i = 0; i < n_citta; i++) {
					buffer_out[i] = path[i]; 
				}

				int send_to = (rank + 1) % size;
				int recv_from = (rank - 1 + size) % size;
				vector<int> buffer_in(n_citta, 0);

				MPI_Request request;
				MPI_Irecv(buffer_in.data(), n_citta, MPI_INT, recv_from, 0, MPI_COMM_WORLD, &request);
				MPI_Send(buffer_out.data(), n_citta, MPI_INT, send_to, 0, MPI_COMM_WORLD);
				MPI_Wait(&request, MPI_STATUS_IGNORE);

				percorso migrato(n_citta, *pop.get_path(0).get_rnd());
				for (int i = 0; i < n_citta; i++) {
					migrato[i] = buffer_in[i];
				}

				int index_replacement = floor(rand.Rannyu(0, n_individui));
				if (migrato.bounds_control()) {
					pop.set_path(index_replacement, migrato);
				} else {
					cerr << "Rank " << rank << ": vincoli non rispettati nel migrato." << endl;
				}
            }
        }
    }

    percorso best = pop.get_best();
    for (int i = 0; i < n_citta; i++) {
        best_path << best[i] << endl;
    }
    best_path.close();
    L1.close();

    // In ogni processo
    vector<int> path_out = pop.get_best().get_path();

    if (rank == 0) {
        vector<vector<int>> ricevuti(size, vector<int>(n_citta));

        for (int k = 0; k < size; k++) { // core 0 riceve da tutti il migliore in una matrice che uso per comparare
            if (k == 0) {
                ricevuti[0] = path_out;
            } else {
                MPI_Recv(ricevuti[k].data(), n_citta, MPI_INT, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        // confronta le lunghezze
        percorso temp(n_citta, *pop.get_path(0).get_rnd());
        for(int j=0; j<n_citta; j++){
            temp[j]=ricevuti[0][j];
        }
        int best_pop = 0;
        double best_length = temp.get_lenght(); // ricostruisci percorso

        for (int k = 1; k < size; k++) {
            for(int j=0; j<n_citta; j++){
            temp[j]=ricevuti[k][j];
            }   
            double len = temp.get_lenght();
            if (len < best_length) {
                best_length = len;
                best_pop = k;
            }
        }

        // salva il migliore
        ofstream best_of_the_best("../OUTPUT/best_of_the_best.dat");
        best_of_the_best << "# Percorso migliore trovato dal core " << best_pop 
                            << " lungo = " << best_length << endl;
        for (int i = 0; i < n_citta; i++) {
            best_of_the_best << ricevuti[best_pop][i] << endl;
        }
        best_of_the_best.close();

    } else {
        // altri processi mandano il percorso
        MPI_Send(path_out.data(), n_citta, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}

/*                                                               
    MPI_Sendrecv(buffer_out.data(), n_citta, MPI_INT, manda_a, 0,
                buffer_in.data(), n_citta, MPI_INT, ricevi_da, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
*/                                                                  
