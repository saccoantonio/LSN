CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio1 : random.cpp random.h funzioni.h funzioni.cpp Esercizio1.cpp 
		$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
		$(CC) -c random.cpp -o random.o $(CFLAGS)
		$(CC) -c Esercizio1.cpp -o Esercizio1.o $(CFLAGS)
		$(CC) funzioni.o random.o Esercizio1.o -o Esercizio1.exe
clean :
		rm *.o *.exe seed.out
