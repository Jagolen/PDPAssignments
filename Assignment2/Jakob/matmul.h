#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int read_input(const char *file_name, double **values);

int write_output(char *file_name, const double *output, int num_values);