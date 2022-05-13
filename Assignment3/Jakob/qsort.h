#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int read_input(const char *file_name, int **values, int size, int *expanded_size);

int write_output(const char *file_name, int *sorted_list, int list_size);

int cmpfunc(const void *a, const void *b);