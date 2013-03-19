// file_handler.c
// provides functions for reading solar system data from a file

#include "file_handler.h"
#include "coordinates.h"
#include <stdlib.h>
#include <stdio.h>

/* Returns a pointer to a solar system specified by the given file.
 * Filename is the path of an ASCII text file formatted as follows:
 * 	system_name number_of_objects coordinate_system
 * 	name mass position velocity
 * 	...
 * 	name mass position velocity */
solar_system *read_solar_system(char *filename) {
	FILE *fp;
	solar_system *s = (solar_system*) malloc (sizeof(solar_system));
	char coordinate_system;
	
	fp = fopen (filename, "r");
	// read in first line of file and check if the coordinate system is specified
	int c_s_specified = fscanf (fp, "%s%i %c\n", s->name, &(s->num_bodies), &coordinate_system);
	// set system_mass to 0 and allocate bodies 
	s->system_mass = 0;
	body *bodies = (body*) malloc (sizeof(body) * s->num_bodies);
	// read in data for each body
	for (int i=0; i<s->num_bodies; i++) {
		fscanf (fp, "%s%lf%lf%lf%lf%lf%lf%lf%lf\n", bodies[i].name, &(bodies[i].mass), &(bodies[i].radius), bodies[i].position,
				bodies[i].position+1, bodies[i].position+2, bodies[i].velocity, bodies[i].velocity+1, bodies[i].velocity+2);
		bodies[i].position_error = 0;	// error assumed to be 0 to start
		bodies[i].velocity_error = 0;
		s->system_mass += bodies[i].mass;
		if (c_s_specified == 3) {	// Checks if the coordinate system has been specified
			if (coordinate_system == 's')	// Checks if the coordinate system is spherical
				to_cartesian (bodies[i].position);	// Converts the position to Cartesian coordinates
		}
	}
	s->bodies = bodies;
	fclose (fp);	// done reading in data
	
	return s;
}
