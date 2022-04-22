#include "matmul.h"

int main(int argc, char **argv){
    
    if (3 != argc) {
        printf("Usage: input_file output_file\n");
        return 1;
	}
    
    //Declaring input data and initial variable
	char *input_name = argv[1];
	char *output_name = argv[2];
    double *timing;
    int rank, size;
    int north, south, east, west, init_start, init_end, matrix_size;
    int period[2] = {1, 1};
    int size_per_dim[2], position[2];
    double **A, **B, **C, **localA, **localB, **localC;

    //Initializing MPI and get the number of processing elements
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    
    /*Creating a cartesian topology and get the ranks. A 2d square topology is assumed
    so processes per dimension is just sqrt(size)*/
    MPI_Comm cart;
    size_per_dim[0] = (int)sqrt(size);
    size_per_dim[1] = (int)sqrt(size);
    MPI_Cart_create(MPI_COMM_WORLD, 2, size_per_dim, period, 1, &cart);

    //Get rank and coordinates
    MPI_Comm_rank(cart, &rank);
    MPI_Cart_coords(cart, rank, 2, position);

    //Find the neighbors
    MPI_Cart_shift(cart, 0, 1, &west, &east);
    MPI_Cart_shift(cart, 1, 1, &north, &south);

    //Process 0 reads the input file
    if(rank == 0){
        	if (0 > (matrix_size = read_input(input_name, &A, &B))) {
		        return 2;
	    }
		timing = (double*)malloc(size*sizeof(double));
    }

    /* initial shift for a and b */
    MPI_Cart_shift(cart, 0, -position[0], &init_start, &init_end);
    MPI_Sendrecv_replace( localA, 1, MPI_INT  init_end, 1,
        init_start, 1 , cart , &status);
    MPI_Cart_shift(cart, 1, -position[1], &init_start, &init_end);
    MPI_Sendrecv_replace( localB, 1, MPI_INT  init_end, 1,
        init_start, 1 , cart , &status);

    /* Compute the local matrix multiplications */
    int local_size = N/size_per_dim[0];
    for(int step; step<local_size;step++){

        for(int i=0;i<local_size;i++)
            for(int k=0;k<local_size;k++)
                for(int j=0;j<local_size;j++)
                    localC[i,j] += localA[i, k]*localB[k, j];
        
        /* the cannon shift, A east, B north */
        MPI_Sendrecv_replace(localA, local_size*local_size, MPI_INT, 
            west, 1, east, 1, cart, &status);
        MPI_Sendrecv_replace(localB, local_size*local_size, MPI_INT, 
            north, 1, south, 1, cart, &status);
    }

    
    
    //MPI_Bcast(&num_values, 1, MPI_INT, 0, cart);







/*     printf("My rank is %d and my coords are (%d, %d). \
Neighbors: North: %d, South: %d, East: %d, West: %d\n",
        rank, position[0], position[1],north, south,east,west);
    int testvector[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    int local[4];
    MPI_Scatter(testvector, 4, MPI_INT, local, 4, MPI_INT, 0, cart);
    printf("Process %d got local values [%d %d %d %d]\n",rank,local[0], 
        local[1], local[2], local[3]); */
    

    MPI_Finalize();

    return 0;
}


int read_input(const char *file_name, int ***A, int ***B) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}
	int matrix_size;
	if (EOF == fscanf(file, "%d", &matrix_size)) {
		perror("Couldn't read matrix size from input file");
		return -1;
	}
	if (NULL == (**A = malloc(matrix_size * sizeof(int)))) {
		perror("Couldn't allocate memory for matrix A");
		return -1;
	}
    for(int i = 0; i<matrix_size; i++){
        if (NULL == (*A[i] = malloc(matrix_size * sizeof(int)))) {
            perror("Couldn't allocate memory for matrix A");
            return -1;
	}

	if (NULL == (**B = malloc(matrix_size * sizeof(int)))) {
		perror("Couldn't allocate memory for matrix B");
		return -1;
	}
    for(int i = 0; i<matrix_size; i++){
        if (NULL == (*B[i] = malloc(matrix_size * sizeof(int)))) {
            perror("Couldn't allocate memory for matrix B");
            return -1;
	    }
    }

	for (int i=0; i<matrix_size; i++) {
		if (EOF == fscanf(file, "%lf", &((*A)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return matrix_size;
}