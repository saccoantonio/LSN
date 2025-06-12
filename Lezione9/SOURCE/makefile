CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio9 : random.cpp random.h funzioni.h funzioni.cpp percorso.h percorso.cpp generazione.h generazione.cpp Esercizio9.cpp 
		$(CC) -c percorso.cpp -o percorso.o $(CFLAGS)
		$(CC) -c generazione.cpp -o generazione.o $(CFLAGS)
		$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
		$(CC) -c random.cpp -o random.o $(CFLAGS)
		$(CC) -c Esercizio9.cpp -o Esercizio9.o $(CFLAGS)
		$(CC) percorso.o generazione.o funzioni.o random.o Esercizio9.o -o Esercizio9.exe

builder : random.cpp random.h funzioni.h funzioni.cpp builder.cpp 
		$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
		$(CC) -c random.cpp -o random.o $(CFLAGS)
		$(CC) -c builder.cpp -o builder.o $(CFLAGS)
		$(CC) funzioni.o random.o builder.o -o builder.exe

clean :
		rm *.o *.exe seed.out