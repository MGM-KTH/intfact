all: main.o
	gcc main.o -lgmp -o main
	rm *.o

main.o: main.c
	gcc -c main.c
