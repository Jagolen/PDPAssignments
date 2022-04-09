#include "stencil.h"

int main(int argc, char **argv) {
	if (5 != argc) {
		printf("Usage: stencil input_file output_file number_of_applications nThreads\n");
		return 1;
	}
	char *input_name = argv[1];
	char *output_name = argv[2];
	int const num_steps = atoi(argv[3]);
	int const nThreads = atoi(argv[4]);

	// Read input file
	double *input;
	int num_values;
	if (0 > (num_values = read_input(input_name, &input))) {
		return 2;
	}

	// Stencil values
	double h = 2.0*PI/num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH/2;
	const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};

	// Start timer
	double start = MPI_Wtime();

	// Allocate data for result
	double *output;
	if (NULL == (output = malloc(num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for output");
		return 2;
	}

	/* Initilize mpi */
	MPI_Status status;
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &size);
	int s=0;

	/* Matrix for local parts of data */
	int const perThread = num_steps/nThreads;
	int **locals = (int **)malloc(nThreads*sizeof(double*));
	for(int i=0;i<nThreads;i++){
		locals[i] = (int *)malloc(perThread*sizeof(double));
	}

	int *b = (int*)malloc(perThread*sizeof(double));

	for(int i=0;i<nThreads;i++)
		for(int j=0;j<perThread;j++){
			locals[i][j] = input[j+i*perThread];
		}

	if(rank==0){
		for(int i=1;i<nThreads;i++){
			MPI_Send( locals[i] , perThread , MPI_DOUBLE , i , i , MPI_COMM_WORLD);
		}
	}
	else{
		MPI_Recv(b , perThread , MPI_DOUBLE , 0 , rank , MPI_COMM_WORLD , &status);
	}

	return 0;
}
