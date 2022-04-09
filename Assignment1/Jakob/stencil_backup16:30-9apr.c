#include "stencil.h"


int main(int argc, char **argv) {
	if (4 != argc) {
		printf("Usage: stencil input_file output_file number_of_applications\n");
		return 1;
	}

    //Declaring input data and initial variable
	char *input_name = argv[1];
	char *output_name = argv[2];
	int num_steps = atoi(argv[3]);
    int size, rank, left, right, num_values;
	const int period = 1;
    double *input;

    //Initializing MPI and defining rank and size
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Creating a cartesian topology and get the ranks and neighbors
    MPI_Comm cart;
    MPI_Cart_create(MPI_COMM_WORLD, 1, &size, &period, 1, &cart);
    MPI_Comm_rank(cart, &rank);

	//Finding the neighbors
	int temp = rank;
	MPI_Cart_shift(cart, 0, -1, &temp, &left);
	MPI_Cart_shift(cart, 0, 1, &temp, &right);
    

    //Rank 0 reads the input data which is then broadcasted to everyone
    
    //TEST CODE REMOVE LATER
    if(rank == 0){
        input = (double*)malloc(12*sizeof(double));
        for(int i = 0; i<12; i++)
            input[i] = i;
        num_values = 12;
    }



/*     if(rank == 0){
        	if (0 > (num_values = read_input(input_name, &input))) {
		        return 2;
	    }
    } */
    MPI_Bcast(&num_values, 1, MPI_INT, 0, cart);

    //This many vals in each process
    int vals_per_pc = num_values/size;

    /*Creating local input and output arrays of same siza as number of values per processor
    with padding of length 2 on both left and right side of the array
    */

    double *localinput = (double*)malloc((4+vals_per_pc)*sizeof(double));
    double *localoutput = (double*)malloc((4+vals_per_pc)*sizeof(double));
    double *l_in = &localinput[2];
    double *l_out = &localoutput[2];

    MPI_Scatter(input, vals_per_pc, MPI_DOUBLE, l_in, vals_per_pc, MPI_DOUBLE, 0, cart);

    for(int i = 0; i<vals_per_pc; i++){
        printf("Element %d in process %d is %f\n",i,rank,l_in[i]);
    }
	printf("For process %d: left = %d, right = %d\n",rank,left,right);

	// Stencil values
	double h = 2.0*PI/num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH/2;
	const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};


    //Free the local arrays

    free(localinput);
    free(localoutput);

    //End of program
    MPI_Finalize();
    return 0;
}


int read_input(const char *file_name, double **values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	if (EOF == fscanf(file, "%d", &num_values)) {
		perror("Couldn't read element count from input file");
		return -1;
	}
	if (NULL == (*values = malloc(num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num_values;
}


int write_output(char *file_name, const double *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%.4f ", output[i])) {
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "\n")) {
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file)) {
		perror("Warning: couldn't close output file");
	}
	return 0;
}
