/**
 * Demo program for MPI_Alltoall
 * Author: Malin Kallen
 */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "print_array.h"

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

	int rank, num_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	int *send_buf = malloc(num_proc*sizeof(int));
	int *recv_buf = calloc(num_proc, sizeof(int));

	for (int i=0; i<num_proc; i++) {
		send_buf[i] = 30 + 10*rank + i;
	}

	printf("proc %d: recv_buf before alltoall: ", rank);
	print_array(recv_buf, num_proc);
	printf("proc %d: send_buf before alltoall: ", rank);
	print_array(send_buf, num_proc);

	MPI_Alltoall(send_buf, 1, MPI_INT, recv_buf, 1, MPI_INT, MPI_COMM_WORLD);

	printf("proc %d: recv_buf after alltoall: ", rank);
	print_array(recv_buf, num_proc);

	free(send_buf);
	free(recv_buf);
	MPI_Finalize();
	return 0;
}
