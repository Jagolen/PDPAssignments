#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/**
*@param x
*@param w
*
*/
void prop(int *x, double *w);
void P_matrix(int *P_vals, int *colms, int *rows);
void row_from_sparse_converter(int r, int *P_vals,
     int *colms, int *rows, int* P);