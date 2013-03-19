// PlanetarySimulation.c
// provides functions for simulating solar systems

#include "PlanetarySimulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#define G 0.0000000000667	// Gravitational constant

// All physics done in Cartesian coordinates
// Units: kg, m, s

static double origin[3] = { 0.0, 0.0, 0.0 };

// Distance between two points
static double dist(double *p1, double *p2) {
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
static body one_step_body(solar_system *s, int index, double step_size) {
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
