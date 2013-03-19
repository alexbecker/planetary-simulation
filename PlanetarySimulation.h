// PlanetarySimulation.h

#ifndef PLANETARYSIMULATION_H
#define PLANETARYSIMULATION_H

typedef struct body {
	char name[20];
	double mass;
	double radius;
	double position[3];
	double velocity[3];
	double position_error;
	double velocity_error;
} body;

typedef struct solar_system {
	char name[20];
	int num_bodies;
	body *bodies;
	double system_mass;
	double system_energy;
} solar_system;

void move_system(solar_system *, double, double, double);

#endif
