// coordinates.c
// provides functions for converting coordinate systems

#include "coordinates.h"
#include <math.h>

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
