/* -----------------------Planetary Simulation----------------------- *
 * This program is intended to compute the positions of astronomical
 * objects given their mass, radius, and initial position and velocity.
 * To do so it numerically solves the n-body problem, keeping track of
 * an upper bound on the error, and shrinks the step size until the
 * error is acceptable. It is called from the command line using:
 * 	./PlanetarySimulation filename endtime max_error [start_step] [c_s]
 * Filename is the path of an ASCII text file formatted as follows:
 * 	system_name number_of_objects coordinate_system
 * 	name mass position velocity
 * 	...
 * 	name mass position velocity
 * Coordinate system is the system of coordinates used to specify the
 * position of each object, and is either "c" for Cartesian or "s" for
 * (mathematical) spherical. The default is Cartesian. Position is the 
 * position of the object in the specified coordinates, seperated by
 * spaces, and velocity is the velocity of the object in Cartesian
 * coordinates. Endtime is the length of time to simulate. Max_error is
 * the largest acceptable error in the final position of any object. The
 * optional argument start_step specifies the step size used initially
 * by the program. Since the error produced by the algorithm is
 * generally much smaller than the upper bound, it may be desirable to
 * specify the step size and set max_error to 1e308, so it is ignored.
 * C_s is the coordinate system which will be used to print the output.
 * ------------------------Alex Becker, 2013------------------------- */

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
