#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int read_input(const char *file_name, int **values);

int write_output(const char *file_name, int *sorted_list, int list_size);