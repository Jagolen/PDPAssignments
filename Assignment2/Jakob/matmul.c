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
    int north, south, east, west, matrix_size;
    int period[2] = {1, 1};
    int size_per_dim[2], position[2];
    int **A, **B, **C, **localA, **localB, **localC;

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
    printf("Hihi\n");
    //Process 0 reads the input file
    if(rank == 0){
        	if (0 > (matrix_size = read_input(input_name, &A, &B))) {
		        return 2;
	    }
		timing = (double*)malloc(size*sizeof(double));
        printf("Matrix A: \n");
        for(int i = 0; i< matrix_size; i++){
            for(int j = 0; j<matrix_size; j++){
                printf("%d ", A[i][j]);
            }
            printf("\n");
        }

        printf("Matrix B: \n");
        for(int i = 0; i< matrix_size; i++){
            for(int j = 0; j<matrix_size; j++){
                printf("%d ", B[i][j]);
            }
            printf("\n");
        }
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
    printf("After matrix size: size = %d\n", matrix_size);


	if (NULL == (*A = malloc(matrix_size * sizeof(int*)))) {
		perror("Couldn't allocate memory for matrix A");
		return -1;
	}
    printf("Bla\n");
    for(int i = 0; i<matrix_size; i++){
        printf("i = %d\n",i);
        if (NULL == (A[i] = malloc(matrix_size * sizeof(int)))) {
            perror("Couldn't allocate memory for matrix A");
            return -1;
        }
	}
    printf("After A\n");
	if (NULL == (*B = malloc(matrix_size * sizeof(int*)))) {
		perror("Couldn't allocate memory for matrix B");
		return -1;
	}
    for(int i = 0; i<matrix_size; i++){
        if (NULL == (B[i] = malloc(matrix_size * sizeof(int)))) {
            perror("Couldn't allocate memory for matrix B");
            return -1;
	    }
    }
    printf("After Mallocs\n");
	for (int i=0; i<matrix_size; i++) {
        for (int j=0; j<matrix_size; j++)
        {
            printf("[%d, %d]\n",i,j);
            if (EOF == fscanf(file, "%d", &((*A)[i][j]))) {
                perror("Couldn't read elements from input file");
                return -1;
            }
        }
	}

    for (int i=0; i<matrix_size; i++) {
        for (int j=0; j<matrix_size; j++)
        {
            if (EOF == fscanf(file, "%d", &((*B)[i][j]))) {
                perror("Couldn't read elements from input file");
                return -1;
            }
        }
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return matrix_size;
}