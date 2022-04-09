/**
 * Demo program for MPI_Scatter
 * Author: Malin Kallen
 */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "print_array.h"

#define ROOT 0

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

	int rank, num_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	int *send_buf = 0;
	if (ROOT == rank) {
		send_buf = malloc(num_proc*sizeof(int));
		for (int i=0; i<num_proc; i++) {
			send_buf[i] = i + 21;
		}
		printf("send_buf on root before scatter: ");
		print_array(send_buf, num_proc);
	}

	int recv_buf = 0;
	printf("proc %d: recv_buf before scatter: %d\n", rank, recv_buf);

	MPI_Scatter(send_buf, 1, MPI_INT, &recv_buf, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	printf("proc %d: recv_buf after scatter: %d\n", rank, recv_buf);

	free(send_buf);
	MPI_Finalize();
	return 0;
}
