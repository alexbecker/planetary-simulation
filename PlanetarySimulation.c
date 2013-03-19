/* -----------------------Planetary Simulation----------------------- *
 * This program is intended to compute the positions of astronomical
 * objects given their mass, radius, and initial position and velocity.
 * To do so it numerically solves the n-body problem, keeping track of
 * an upper bound on the error, and shrinks the step size until the
 * error is acceptable. It is called from the command line using:
 * 	./PlanetarySimulation filename endtime max_error [start_step]
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
 * ------------------------Alex Becker, 2013------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#define G 0.0000000000667	// Gravitational constant

/*---- PHYSICS SIMULATOR ----*/
// All physics done in Cartesian coordinates
// Units: kg, m, s

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

static double origin[3] = { 0.0, 0.0, 0.0 };

// Distance between two points
double dist(double *p1, double *p2) {
	return sqrt (pow (p1[0] - p2[0], 2) + pow (p1[1] - p2[1], 2) + pow (p1[2] - p2[2], 2));
}

// Moves a body one step
// The approximation assumes that acceleration is constant
// Keeps track of the error induced
/* The true error in velocity EACH STEP is bounded above by
 *  step_size*|acceleration_error|+1/2*step_size^2*M 
 * where M is an upper bound on the magnitude of the jerk.
 * The true error in position EACH STEP is bounded above by
 *  step_size*|velocity_error|+1/6*step_size^3*M
 * This is the basis for computing error. */
body one_step_body(solar_system *s, int index, double step_size) {
	body p = s->bodies[index];
	double cur_dist, min_dist = DBL_MAX, max_velocity = 0, t, new_min_dist, falling_speed, M, multiplier, new_velocity[3];
	for (int i=0; i<s->num_bodies; i++) {
		if (i != index) {
			cur_dist = dist (p.position, s->bodies[i].position);	// Distance between bodies
			if (cur_dist <= p.radius + s->bodies[i].radius) {		// Checks for collisions
				printf ("COLLISION DETECTED BETWEEN %s AND %s\n", p.name, s->bodies[i].name);
				exit (0);					// Program cannot simulate collisions, so terminates
			}
			// Lower bound on distance between bodies, taking into account existing errors in position
			min_dist = fmin (min_dist, cur_dist - p.position_error - s->bodies[i].position_error);	
		}
		// Upper bound on speed of any object in the system, taking into account existing errors in velocity
		max_velocity = fmax (max_velocity, dist (origin, s->bodies[i].velocity) + s->bodies[i].velocity_error);
	}
	if (min_dist <= 0) {			// Makes following assignments safe
		p.position_error = DBL_MAX / 2;		// Divide by two to avoid overflow
	} else {	// The following code exists solely to calculate M
		/* t is a lower bound for the amount of time it would take for the system to collide if all of its mass were concentrated in the two
		 * closest bodies, and they were both travelling towards each other at the maximum velocity of any body in the system. */
		t = min_dist / (2 * sqrt (pow (s->system_mass, 2) / (p.mass * (s->system_mass - p.mass)) * max_velocity + 4 * G * s->system_mass / min_dist));
		if (step_size >= t)	// t is useless if it is less than step_size
			p.position_error = DBL_MAX / 2;
		// new_min_dist is the distance between two such objects after step_size
		new_min_dist = (1 - step_size / t) * min_dist;
		// falling_speed is the rate at which the distance between them would be decreasing
		falling_speed = sqrt (2 * G * (1 / new_min_dist - 1 / min_dist) / p.mass + pow (2 * max_velocity, 2));
		M = s->system_mass * G * falling_speed / pow (new_min_dist, 3);	// Upper bound on jerk
	}
	for (int i=0; i<3; i++)
		new_velocity[i] = p.velocity[i];
	for (int j=0; j<s->num_bodies; j++) {
		if (j != index) {
			cur_dist = dist (p.position, s->bodies[j].position); // Re-used this variable cuz why not
			multiplier = step_size * s->bodies[j].mass * G / pow (cur_dist, 3);
			if (cur_dist <= p.position_error + s->bodies[j].position_error)	// makes next assignment accurate
				p.position_error = DBL_MAX / 2;					// Divide by two to avoid overflow
			// Next line adds step_size*|acceleration_error|
			p.velocity_error += step_size * s->bodies[j].mass * G * (1 / pow (cur_dist - p.position_error - s->bodies[j].position_error, 2) - 1 / pow (cur_dist, 2));
			for (int i=0; i<3; i++)
				new_velocity[i] += (s->bodies[j].position[i] - p.position[i]) * multiplier;
		}
	}
	for (int i=0; i<3; i++) {
		p.position[i] += step_size * (p.velocity[i] + new_velocity[i]) / 2;	// Trapezoidal rule
		p.velocity[i] = new_velocity[i];
	}
	p.velocity_error += pow (step_size, 2) * M / 2;
	p.position_error += step_size * p.velocity_error + pow (step_size, 3) * M / 6;
	return p;
}

// Moves the entire solar_system for specified amount of time
// Guaranteed to be accurate within the the given max_error
void move_system(solar_system *s, double max_error, double endtime, double step_size) {
	double time, error;
	body *old_bodies, *new_bodies;
	old_bodies = malloc (sizeof (body) * s->num_bodies);
	for (int i=0; i<s->num_bodies; i++)
		old_bodies[i] = s->bodies[i];
	while (1) {		// Loops until error < max_error
		time = 0;
		error = 0;
		while (time < endtime) {
			new_bodies = malloc (sizeof (body) * s->num_bodies);
			for (int i=0; i<s->num_bodies; i++)
				new_bodies[i] = one_step_body (s, i, fmin (endtime - time, step_size));	// Step not more than remaining time
			free (s->bodies);
			s->bodies = new_bodies;
			time += step_size;
		}
		for (int i=0; i<s->num_bodies; i++)
			error = fmax (error, s->bodies[i].position_error);
		printf ("Step = %fs, Error < %fm\n", step_size, error);	// So the user knows what's happening
		if (error < max_error)
			break;				// Break if error < max_error
		for (int i=0; i<s->num_bodies; i++)
			s->bodies[i] = old_bodies[i];	// Otherwise reset bodies and try again
		if (error >= DBL_MAX / 2) {		// Catches huge errors
			step_size /= 10;
		} else {
			step_size *= .95 * max_error / error;	// Change step size, relying on approx linearity of error
		}
	}
	return;
}

/*---- I/O PROCESSING ----*/

// Converts Cartesian coordinates to spherical coordinates
// Spherical coordinates use mathematical notation
void to_spherical(double *cartesian) {
	double spherical[3];
	spherical[0] = sqrt (pow (cartesian[0], 2) + pow (cartesian[1], 2) + pow (cartesian[2], 2));
	spherical[1] = atan2 (cartesian[1], cartesian[0]);
	spherical[2] = atan2 (cartesian[2], sqrt (pow (cartesian[0], 2) + pow (cartesian[1], 2))) + 1.5707963267;
	for (int i=0; i<3; i++)
		cartesian[i] = spherical[i];
	return;
}

// Converts spherical coordinates to Cartesian coordinates
void to_cartesian(double *spherical) {
	double cartesian[3];
	cartesian[0] = sin (spherical[2]) * cos (spherical[1]) * spherical[0];
	cartesian[1] = sin (spherical[2]) * sin (spherical[1]) * spherical[0];
	cartesian[2] = cos (spherical[2]) * spherical[0];
	for (int i=0; i<3; i++)
		spherical[i] = cartesian[i];
	return;
}

int main(int argc, char **argv) {
	FILE *fp;
	double endtime, max_error, start_step;
	body *bodies;
	solar_system *s;
	char coordinate_system;
	int c_s_specified;
	s = (solar_system*) malloc (sizeof(solar_system));
	fp = fopen (argv[1], "r");
	endtime = strtod (argv[2], (char **)0);
	max_error = strtod (argv[3], (char **)0);
	if (argc > 4) {						// Checks if start_step specified
		start_step = strtod (argv[4], (char **)0);	// Reads specified start_step
	} else {
		start_step = fmin (endtime / 1000, 100);	// Otherwise uses default value
	}
	c_s_specified = fscanf (fp, "%s%i %c\n", s->name, &(s->num_bodies), &coordinate_system);
	s->system_mass = 0;
	bodies = (body*) malloc (sizeof(body) * s->num_bodies);
	for (int i=0; i<s->num_bodies; i++) {
		fscanf (fp, "%s%lf%lf%lf%lf%lf%lf%lf%lf\n", bodies[i].name, &(bodies[i].mass), &(bodies[i].radius), bodies[i].position,
				bodies[i].position+1, bodies[i].position+2, bodies[i].velocity, bodies[i].velocity+1, bodies[i].velocity+2);
		bodies[i].position_error = 0;
		bodies[i].velocity_error = 0;
		s->system_mass += bodies[i].mass;
		if (c_s_specified == 3) {	// Checks if the coordinate system has been specified
			if (coordinate_system == 's')	// Checks if the coordinate system is spherical
				to_cartesian (bodies[i].position);	// Converts the position to Cartesian coordinates
		}
	}
	fclose (fp);
	s->bodies = bodies;
	printf ("Simulating %s:\n", s->name);
	move_system (s, max_error, endtime, start_step);	// Performs simulation
	printf ("No collisions have occured.\n");	// If collisions have occured, program already terminated
	for (int i=0; i<s->num_bodies; i++) {
		if (c_s_specified == 3) {
			if (coordinate_system == 's') {
				to_spherical (s->bodies[i].position); // Converts the position back to spherical coordinates
				printf ("%s is %fm from the center, at azimuthal angle %f and polar angle %f, with velocity vector (%fm/s, %fm/s, %fm/s).\n",
						s->bodies[i].name, s->bodies[i].position[0], s->bodies[i].position[1], s->bodies[i].position[2], 
						s->bodies[i].velocity[0], s->bodies[i].velocity[1], s->bodies[i].velocity[2]);
			} else {
				printf ("%s is located at (%fm, %fm, %fm), with velocity vector (%fm/s, %fm/s, %fm/s).\n",
						s->bodies[i].name, s->bodies[i].position[0], s->bodies[i].position[1], s->bodies[i].position[2], 
						s->bodies[i].velocity[0], s->bodies[i].velocity[1], s->bodies[i].velocity[2]);
			}
		} else {
			printf ("%s is located at (%fm, %fm, %fm), with velocity vector (%fm/s, %fm/s, %fm/s).\n",
					s->bodies[i].name, s->bodies[i].position[0], s->bodies[i].position[1], s->bodies[i].position[2], 
					s->bodies[i].velocity[0], s->bodies[i].velocity[1], s->bodies[i].velocity[2]);
		}
	}
	return 0;
}
