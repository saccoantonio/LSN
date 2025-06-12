CXX = mpic++
CXXFLAGS = -Wall -O3 -std=c++11

# File oggetto
OBJ = percorso.o generazione.o funzioni.o random.o Esercizio10.o

# Target principale
Esercizio10: $(OBJ)
	$(CXX) $(OBJ) -o Esercizio10.exe

# Regole per i singoli oggetti
percorso.o: percorso.cpp percorso.h
	$(CXX) -c percorso.cpp -o percorso.o $(CXXFLAGS)

generazione.o: generazione.cpp generazione.h
	$(CXX) -c generazione.cpp -o generazione.o $(CXXFLAGS)

funzioni.o: funzioni.cpp funzioni.h
	$(CXX) -c funzioni.cpp -o funzioni.o $(CXXFLAGS)

random.o: random.cpp random.h
	$(CXX) -c random.cpp -o random.o $(CXXFLAGS)

Esercizio10.o: Esercizio10.cpp random.h funzioni.h percorso.h generazione.h
	$(CXX) -c Esercizio10.cpp -o Esercizio10.o $(CXXFLAGS)

# Pulizia
clean:
	rm -f *.o *.exe seed.out
