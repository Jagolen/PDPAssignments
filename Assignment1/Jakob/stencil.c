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
    int size, rank, west, east, num_values;
	const int period = 1;
    double *in_out;

    //Initializing MPI and defining rank and size
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Creating a cartesian topology and get the ranks and neighbors
    MPI_Comm cart;
    MPI_Cart_create(MPI_COMM_WORLD, 1, &size, &period, 1, &cart);
    MPI_Comm_rank(cart, &rank);

	//Finding the neighbors
	MPI_Cart_shift(cart, 0, 1, &west, &east);

	//Statuses for left and right
	MPI_Status w_stat, e_stat;
    
    //Rank 0 reads the input data which is then broadcasted to everyone
    
    //TEST CODE REMOVE LATER
/*     if(rank == 0){
        input = (double*)malloc(12*sizeof(double));
        for(int i = 0; i<12; i++)
            input[i] = i;
        num_values = 12;
    } */



    if(rank == 0){
        	if (0 > (num_values = read_input(input_name, &in_out))) {
		        return 2;
	    }
    }
    MPI_Bcast(&num_values, 1, MPI_INT, 0, cart);

    //This many vals in each process
    int vals_per_pc = num_values/size;

    /*Creating local input and output arrays of same size as number of values per processor
    with padding of length 2 on both left and right side of the array
    */

    double *buflocalinput = (double*)malloc((4+vals_per_pc)*sizeof(double));
    double *l_out = (double*)malloc((vals_per_pc)*sizeof(double));
    double *l_in = &buflocalinput[2];

    MPI_Scatter(in_out, vals_per_pc, MPI_DOUBLE, l_in, vals_per_pc, MPI_DOUBLE, 0, cart);

	// Stencil values
	double h = 2.0*PI/num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH/2;
	const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};


	//Creating the data type to send two elements to a neighbor
	MPI_Datatype edge;
	MPI_Type_contiguous(EXTENT, MPI_DOUBLE, &edge);
	MPI_Type_commit(&edge);

	//Setting up request handles
	MPI_Request w_send_req, e_send_req, w_rec_req, e_rec_req;
	
	//Sending and recieving data between neighbors
	MPI_Isend(&l_in[0],1,edge,west,west,cart,&w_send_req);
	MPI_Isend(&l_in[vals_per_pc-EXTENT],1,edge,east,east,cart,&e_send_req);
	MPI_Irecv(&buflocalinput[0],1,edge,west,rank,cart,&w_rec_req);
	MPI_Irecv(&buflocalinput[vals_per_pc+EXTENT],1,edge,east,rank,cart,&e_rec_req);

	//Performing the stencil on the elements in the middle first

	for(int i = EXTENT; i<vals_per_pc-EXTENT; i++){
		double result = 0;
		for(int j = 0; j<STENCIL_WIDTH;j++){
			int index = i-EXTENT+j;
			result += STENCIL[j]*l_in[index];
		}
		l_out[i] = result;
	}

	//Waiting for the data from the west to be recieved and then perform the stencil on the left elements
	MPI_Wait(&w_rec_req, &w_stat);
	
	for(int i = 0; i<EXTENT; i++){
		double result = 0;
		for(int j = 0; j<STENCIL_WIDTH;j++){
			int index = i-EXTENT+j;
			result += STENCIL[j]*l_in[index];
		}
		l_out[i] = result;
	}

	//Waiting for the data from the east to be recieved and then perform the stencil on the right elements
	MPI_Wait(&e_rec_req, &e_stat);
	for(int i = vals_per_pc-EXTENT; i<vals_per_pc; i++){
		double result = 0;
		for(int j = 0; j<STENCIL_WIDTH;j++){
			int index = i-EXTENT+j;
			result += STENCIL[j]*l_in[index];
		}
		l_out[i] = result;
	}

	if(rank == 0) printf("process 0 got %f and %f from %d and %f and %f from %d\n",
		buflocalinput[0], buflocalinput[1], west, buflocalinput[vals_per_pc], 
			buflocalinput[vals_per_pc+1], east);
	//Copying values from output to input so the stencil can be applied again
	for(int i = 0; i<vals_per_pc; i++){
		l_in[i] = l_out[i];
	}


	//Gathering each contribution from the processess
	MPI_Gather(l_out, vals_per_pc, MPI_DOUBLE, in_out, vals_per_pc, MPI_DOUBLE, 0, cart);

	//Process 0 creates an output file
	if (0 != write_output(output_name, in_out, num_values)) {
		return 2;
	}

	//printf("Process %d recieved %f and %f from process %d\n",rank,buflocalinput[0], buflocalinput[1],west);
	//printf("Process %d recieved %f and %f from process %d\n",rank,buflocalinput[vals_per_pc+EXTENT], buflocalinput[vals_per_pc+EXTENT+1],east);
    
	
	
	//Free the local arrays
    free(buflocalinput);
    free(l_out);
	if(rank == 0) free(in_out);
    
	//End of program
	MPI_Type_free(&edge);
	MPI_Comm_free(&cart);
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
