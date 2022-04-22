#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int read_input(const char *file_name, int ***A, int ***B);

int write_output(char *file_name, const double *output, int num_values);