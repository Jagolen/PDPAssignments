#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/**
*@param x
*@param w
*@param P_vals
*@param file_name
*@param results
*@param N
*/
void prop(int *x, double *w);
void P_matrix(int *P_vals);
int write_output(const char *file_name, int *results, int N);

