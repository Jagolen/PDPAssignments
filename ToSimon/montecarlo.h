#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/**
*@param x State vector Should be of length 7!
*@param w Result vector (propensities). Should be of length 15!
*/
void prop(int *x, double *w);

int write_output(const char *file_name, int *result, int num_bins, int *bins);