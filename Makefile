# Makefile for Planetary Simulation

CC = gcc
CFLAGS = -std=gnu99
LOADLIBES = -lm
SRC = PlanetarySimulation.c file_handler.c coordinates.c launcher.c
OBJ = PlanetarySimulation.o file_handler.o coordinates.o launcher.o

launcher: ${OBJ}

file_handler: file_handler.o coordinates.o

coordinates: coordinates.o

PlanetarySimulation: PlanetarySimulation.o

clean:
	rm -f ${OBJ} launcher
	
depend:
	./remove_depends.sh
	gcc -MM ${SRC} >> Makefile
	
# dependencies listed by gcc -MM
PlanetarySimulation.o: PlanetarySimulation.c PlanetarySimulation.h
file_handler.o: file_handler.c file_handler.h coordinates.h \
 PlanetarySimulation.h
coordinates.o: coordinates.c coordinates.h
launcher.o: launcher.c PlanetarySimulation.h file_handler.h coordinates.h
