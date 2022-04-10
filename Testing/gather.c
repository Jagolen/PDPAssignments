/**
 * Demo program for MPI_Gather
 * Author: Malin Kallen
 */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include <print_array.h>

#define ROOT 0

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

	int rank, num_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	int send_buf = rank + 11;
	printf("proc %d: send_buf before gather: %d\n", rank, send_buf);

	int *recv_buf = 0;
	if (ROOT == rank) {
		recv_buf = calloc(num_proc, sizeof(int));
		printf("recv_buf on root before gather: ");
		print_array(recv_buf, num_proc);
	}

	MPI_Gather(&send_buf, 1, MPI_INT, recv_buf, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	if (ROOT == rank) {
		printf("recv_buf on root after gather: ");
		print_array(recv_buf, num_proc);
	}

	free(recv_buf);
	MPI_Finalize();
	return 0;
}
