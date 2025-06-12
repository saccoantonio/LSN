CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio3 : random.cpp random.h funzioni.h funzioni.cpp Esercizio3.cpp 
		$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
		$(CC) -c random.cpp -o random.o $(CFLAGS)
		$(CC) -c Esercizio3.cpp -o Esercizio3.o $(CFLAGS)
		$(CC) funzioni.o random.o Esercizio3.o -o Esercizio3.exe

clean :
		rm *.o *.exe seed.out
