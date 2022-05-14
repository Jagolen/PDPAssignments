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
	MPI_Status status;

	//Rank 0 reads the input list and broadcast number of elements
    if(rank == 0){
        if (0 > (list_size = read_input(input_name, &values))) return 2;
		timing = (double*)malloc(size*sizeof(double));
    }
	MPI_Bcast(&list_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Finding the number of elements per local array and the remainder if it isn't divisible.
	int nr_local_values = list_size/size;
	int extras = list_size%size;

	//Setting up values for a Scatterv
	int displacement[size];
	int count_send[size];
	int displacement_now = 0;

	for(int i = 0; i<size; i++){
		displacement[i] = displacement_now;
		if(i == 0){
			count_send[i] = nr_local_values + extras;
			displacement_now += (nr_local_values + extras);
		}
		else{
			count_send[i] = nr_local_values;
			displacement_now += nr_local_values;
		}
	}

	

	//The extra elements are put into rank 0
	if(rank == 0) nr_local_values += extras;

	printf("%d LocalN = %d\n", rank, nr_local_values);

	//Creating local lists
	if(nr_local_values == 0) local_values = (int*)malloc(sizeof(int));
	else local_values = (int*)malloc(nr_local_values*sizeof(int));

	//To merge the lists
	merge_array = (int*)malloc(sizeof(int));

	//To recieve elements from other processes
	recieve_array = (int*)malloc(sizeof(int));


	//The values are scattered to the local lists and a timer starts
	
    MPI_Scatterv(values, count_send, displacement, MPI_INT, local_values, 
		nr_local_values, MPI_INT, 0, hyper);

	//MPI_Scatter(values, nr_local_values, MPI_INT, local_values, nr_local_values, MPI_INT, 0, hyper);
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

		//Get the size of the sub group and the new ranks
		MPI_Comm_rank(hyper_sub, &sub_rank);
		MPI_Comm_size(hyper_sub, &sub_size);

		//We get the median for every process
		int median, chosen_median;

		if(nr_local_values == 0) median = 0;
		else median = local_values[nr_local_values/2];

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
				}
				exit(-1);
			}
		}

		//Next we find which element in the local array is the pivot element
		int pivot_element = 0;
		while(pivot_element < nr_local_values && local_values[pivot_element]<chosen_median) pivot_element++;

		/*
		The lowest in the lower half sends to the lowest in the upper half, the second
		lowest in the lower half to the second lowest in the upper half and so on.
		*/

		int partner = (sub_rank + (sub_size/2))%sub_size;
		printf("Partner of %d is %d\n", rank, partner);
		int sendcount, recievecount, nr_merge_elements, mergealloc;

		//The lower half
		if(sub_rank < (sub_size/2)){
			/*
			The lower half sends the elements above the pivot, the partner sends how many elements
			the current process will get and the current process sends how many the partner will get.
			*/
			sendcount = nr_local_values - pivot_element;
			MPI_Sendrecv(&sendcount, 1, MPI_INT, partner, 1, &recievecount, 1, MPI_INT, partner, 1, hyper_sub, &status);

			//Reallocates as much memory as needed in the recieve and merge arrays
			if(recievecount == 0) recieve_array = (int*)realloc(recieve_array, sizeof(int));
			else recieve_array = (int*)realloc(recieve_array, recievecount*sizeof(int));
			nr_merge_elements = pivot_element+recievecount;
			mergealloc = nr_merge_elements*sizeof(int);
			if(mergealloc == 0) merge_array = (int*)realloc(merge_array, sizeof(int));
			else merge_array = (int*)realloc(merge_array, mergealloc);

			//printf("%d gets %d from %d and sends %d\n", rank, recievecount, partner, sendcount);
			//Send and recieve the data

			//printf("%d sending\n", rank);
			MPI_Sendrecv(local_values+pivot_element, nr_local_values-pivot_element, 
				MPI_INT, partner, 1, recieve_array, recievecount, MPI_INT, partner, 1, 
				hyper_sub, &status);
			//printf("%d recieved\n", rank);
		}

		//The upper half
		else{
			//The upper half sends the elements up to the pivot
			sendcount = pivot_element;
			MPI_Sendrecv(&sendcount, 1, MPI_INT, partner, 1, &recievecount, 1, MPI_INT, partner, 1, hyper_sub, &status);

			//Again memory for the arrays are reallocated
			if(recievecount == 0) recieve_array = (int*)realloc(recieve_array, sizeof(int));
			else recieve_array = (int*)realloc(recieve_array, recievecount*sizeof(int));
			nr_merge_elements = nr_local_values-pivot_element+recievecount;
			mergealloc = nr_merge_elements*sizeof(int);
			if(mergealloc == 0) merge_array = (int*)realloc(merge_array, sizeof(int));
			else merge_array = (int*)realloc(merge_array, mergealloc);

			//Send and recieve the data
			MPI_Sendrecv(local_values, pivot_element, MPI_INT, partner, 1, 
			recieve_array, recievecount, MPI_INT, partner, 1, hyper_sub, &status);
		}

		//Now the local values are merged with the recieved values and sorted
	
		//Initializing parameters
		int val, val_max;
		int rec = 0;
		int rec_max = recievecount;
		int merge = 0;

		//Lower half, local values are checked from 0 to the pivot element
		if(sub_rank < (sub_size/2)){
			val = 0;
			val_max = pivot_element;
		}

		//Else local values are checked from the pivot element to the end of the list
		else{
			val = pivot_element;
			val_max = nr_local_values;
		}

		/*
		The values of the recieve array and the local array are compared, whichever is lowest
		is inserted into the merge array and the array the value was taken from
		advances one element. This ensures that the list will be sorted.
		*/

		while(val<val_max && rec<rec_max){
			if(local_values[val] < recieve_array[rec]){
				merge_array[merge] = local_values[val];
				merge++;
				val++;
			}
			else{
				merge_array[merge] = recieve_array[rec];
				merge++;
				rec++;
			}
		}

		//If one array has finished, the data from the other array is copied
		while(val<val_max){
			merge_array[merge] = local_values[val];
			merge++;
			val++;
		}
		while(rec<rec_max){
			merge_array[merge] = recieve_array[rec];
			merge++;
			rec++;
		}

		/*
		With all the data sorted in the merge array, the required memory is allocated
		for the local array and the data is copied from the merge array
		*/

		if(mergealloc == 0) local_values = (int*)realloc(local_values,sizeof(int));
		else local_values = (int*)realloc(local_values,mergealloc);

		//The new number of elements in the local array is mergealloc
		nr_local_values = nr_merge_elements;

		for(int i = 0; i < nr_local_values; i++){
			local_values[i] = merge_array[i];
		}
		printf("%d NEW N = %d\n", rank, nr_local_values);
		//Free the split communicator so a new one can be used in the next iteration
		MPI_Comm_free(&hyper_sub);
	}

	for(int i = 0; i<size; i++){
		if(rank == i){
			printf("Rank %d has: [", rank);
			for(int j = 0; j<nr_local_values;j++){
				printf("%d ", local_values[j]);
			}
			printf("]\n");
			MPI_Barrier(hyper);
		}
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
int read_input(const char *file_name, int **values) {
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

	if (NULL == (*values = malloc(num_values * sizeof(int)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}

	for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%d", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
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