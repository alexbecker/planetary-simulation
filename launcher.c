// launcher.c
// launches the PlanetarySimulation program

#include "PlanetarySimulation.h"
#include "file_handler.h"
#include "coordinates.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {
	double endtime, max_error, start_step;
	char output_c_s = 'c';	// defaults to cartesian
	// process command line args
	solar_system *s = read_solar_system(argv[1]);
	endtime = strtod (argv[2], (char **)0);
	max_error = strtod (argv[3], (char **)0);
	if (argc > 4) {
		if (sscanf(argv[4], "%f", &start_step)) {
			start_step = strtod (argv[4], (char **)0);	// Reads specified start_step
		} else {
			start_step = fmin (endtime / 1000, 100);	// Otherwise uses default value for start_step
			output_c_s = argv[4][0];					// 4th argument must be the coordinate system
		}
		if (argc > 5) {
			output_c_s = argv[5][0];
		}
	} else {
		start_step = fmin (endtime / 1000, 100);
	}
	// perform simulation
	printf ("Simulating %s:\n", s->name);
	move_system (s, max_error, endtime, start_step);
	printf ("No collisions have occured.\n");	// If collisions have occured, program already terminated
	for (int i=0; i<s->num_bodies; i++) {
		if (output_c_s == 's') {
			to_spherical (s->bodies[i].position); // Converts the position to spherical coordinates
			printf ("%s is %fm from the center, at azimuthal angle %f and polar angle %f, with velocity vector (%fm/s, %fm/s, %fm/s).\n",
						s->bodies[i].name, s->bodies[i].position[0], s->bodies[i].position[1], s->bodies[i].position[2], 
						s->bodies[i].velocity[0], s->bodies[i].velocity[1], s->bodies[i].velocity[2]);
		} else {
			printf ("%s is located at (%fm, %fm, %fm), with velocity vector (%fm/s, %fm/s, %fm/s).\n",
					s->bodies[i].name, s->bodies[i].position[0], s->bodies[i].position[1], s->bodies[i].position[2], 
					s->bodies[i].velocity[0], s->bodies[i].velocity[1], s->bodies[i].velocity[2]);
		}
	}
}
