CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio2 : random.cpp random.h funzioni.h funzioni.cpp Esercizio2.cpp 
		$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
		$(CC) -c random.cpp -o random.o $(CFLAGS)
		$(CC) -c Esercizio2.cpp -o Esercizio2.o $(CFLAGS)
		$(CC) funzioni.o random.o Esercizio2.o -o Esercizio2.exe
clean :
		rm *.o *.exe seed.out
