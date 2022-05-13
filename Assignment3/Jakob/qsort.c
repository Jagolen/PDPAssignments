#include "qsort.h"
int main(int argc, char **argv){
    
    if (4 != argc) {
        printf("Usage: input_file output_file method\n");
        return 1;
	}

	//Parameters
	char *input_name = argv[1];
	char *output_name = argv[2];
    double *timing;
    int rank, size, list_size, expanded_size;
    int *values, *local_values, *merge_array, *recieve_array;
	int method = atoi(argv[3]);

	//Initializing MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//Determining dimensions and parameters for cart_create
	int dims = 0;
	int elements = 1;
	while(elements != size){
		elements *= 2;
		dims++;
	}
	int periodic[dims];
	int proc_per_dim[dims];
	for(int i = 0; i<dims; i++){
		periodic[i] = 0;
		proc_per_dim[i] = 2;
	}

	//Creating a hypercube communicator
	MPI_Comm hyper;
	MPI_Cart_create(MPI_COMM_WORLD, dims, proc_per_dim, periodic, 1, &hyper);
	MPI_Comm_rank(hyper, &rank);

	//Rank 0 reads the input list and broadcast number of elements
    if(rank == 0){
        if (0 > (list_size = read_input(input_name, &values, size, &expanded_size))) return 2;
		timing = (double*)malloc(size*sizeof(double));
		printf("Size = %d, Expanded = %d\n",list_size, expanded_size);
    }
	MPI_Bcast(&list_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&expanded_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Creating local lists
	local_values = (int*)malloc(expanded_size*sizeof(int));

	//To merge the lists
	merge_array = (int*)malloc(expanded_size*sizeof(int));

	//To recieve elements from other processes
	recieve_array = (int*)malloc(expanded_size*sizeof(int));

	//number of elements in the local list
	int nr_local_values = expanded_size/size;

	//The values are scattered to the local lists and a timer starts
	MPI_Scatter(values, nr_local_values, MPI_INT, local_values, nr_local_values, MPI_INT, 0, hyper);
	MPI_Wtime();

	//The first step is to sort each local array with the built in qsort function
	qsort(local_values, nr_local_values, sizeof(int), cmpfunc);

	//We do log_2(p) iterations = log_10(p)/log_10(2)
	int number_of_iterations = round(log(size)/log(2));

	//Now we start iterating
	for(int iteration = 0; iteration<number_of_iterations; iteration++){
		/*Based on the iteration, we split the communicator into groups,
		a single group for iteration 0, 2 groups for iteration 1, 4 groups for iteration 2 and so on
		*/
		MPI_Comm hyper_sub;
		int sub_rank, sub_size;
		int group = pow(2,iteration)*rank/size;
		MPI_Comm_split(hyper, group, rank, &hyper_sub);
		MPI_Status status;

		//Get the size of the sub group and the new ranks
		MPI_Comm_rank(hyper_sub, &sub_rank);
		MPI_Comm_size(hyper_sub, &sub_size);

		//We get the median for every process
		int median, chosen_median;
		median = local_values[nr_local_values/2];

		/*
		The different methods. Method 1 uses the median from process 0, Method 2 takes
		the median of all medians and method 3 uses the mean value of all medians.
		*/
		if(method == 1){
			chosen_median = median;
			MPI_Bcast(&chosen_median, 1, MPI_INT, 0, hyper_sub);
			if(sub_rank == 0){
				printf("Chosen median = %d\n", chosen_median);
			}
		}

		//If we use method 2 or 3 we have to create a list of all medians
		else{
			int medians[sub_size];

			//All medians are gathered in a list
			MPI_Gather(&median, 1, MPI_INT, medians, 1, MPI_INT, 0, hyper_sub);
			if(method == 2){
				if(sub_rank == 0){
					// For method 2 we sort the array of medians and pick the median of that
					qsort(medians, sub_size, sizeof(int), cmpfunc);
					chosen_median = medians[sub_size/2];
					MPI_Bcast(&chosen_median, 1, MPI_INT, 0, hyper_sub);
					printf("Chosen median = %d\n", chosen_median);
				}
			}
			else if(method == 3){
				if(sub_rank == 0){
					//For method 3 we take the mean value of the medians, no sorting required
					double temp_median = 0;
					for(int i = 0; i<sub_size; i++){
						temp_median += (medians[i]/(double)sub_size);
						printf("Median = %d ", medians[i]);
					}
					chosen_median = (int)temp_median;
					MPI_Bcast(&chosen_median, 1, MPI_INT, 0, hyper_sub);
					printf("Chosen median = %d\n", chosen_median);
				}
			}

			//If an invalid method is used, just give an error and exit the program
			else{
				if(rank == 0){
					printf("Invalid method\n");
					exit(-1);
				}
			}
		}

		//Next we find which element in the local array


		//Free the split communicator so a new one can be used in the next iteration
		MPI_Comm_free(&hyper_sub);
	}


    /* writing the output */
    int output_status = 1; // 1 if active output. 0 inactive
    if(rank == 0 && output_status==1){
		if (0 != write_output(output_name, values, list_size)) {
			return 2;
		}
    }


	MPI_Finalize();
    return 0;
}


//Function for reading the input
int read_input(const char *file_name, int **values, int size, int *expanded_size) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	int expanded_values;
	if (EOF == fscanf(file, "%d", &num_values)) {
		perror("Couldn't read element count from input file");
		return -1;
	}

	int remainder = size - (num_values%size);
	if(num_values%size == 0) remainder = 0;
	expanded_values = num_values + remainder;
	printf("Expanded = %d\n", expanded_values);

	if (NULL == (*values = malloc(expanded_values * sizeof(int)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for(int i = 0; i< remainder; i++){
		(*values)[i] = 0;
	}

	for (int i=remainder; i<expanded_values; i++) {
		if (EOF == fscanf(file, "%d", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	*expanded_size = expanded_values;
	return num_values;
}


//Function for writing the output
int write_output(const char *file_name, int *sorted_list, int list_size) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < list_size; i++) {
		if (0 > fprintf(file, "%d ", sorted_list[i])) {
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

//The standard comparison function for qsort
int cmpfunc(const void *a, const void *b){
	return (*(int*)a - *(int*)b);
}