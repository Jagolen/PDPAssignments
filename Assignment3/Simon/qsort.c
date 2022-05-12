#include "qsort.h"


int cmpfunc (const void * a, const void * b) {
	return ( *(int*)a - *(int*)b );
}

int main(int argc, char **argv){
    
    if (4 != argc) {
        printf("Usage: input_file output_file, pivot type 1, 2 or 3\n");
        return 1;
	}
	char *input_name = argv[1];
	char *output_name = argv[2];
	const int pivot = atoi(argv[3]);
    double *timing;
    int rank, size, list_size, extras, local_size;
    int *values, *out, *local_values;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    MPI_Comm hyperCube;
    int nDim = (int)log2(size);
    printf("Dimentions = %d\n",nDim);
    int local_sizes[nDim], position[nDim], period[nDim];
    for(int i=0;i<nDim;i++) {
        local_sizes[i]=2;
        period[i]=1;
    }
    MPI_Cart_create(MPI_COMM_WORLD, nDim, local_sizes, period, 1, &hyperCube);
    int my_coords[nDim];
    MPI_Cart_coords(hyperCube, rank, nDim, my_coords);
    
	/* reading input where values is the list and list_size is the 
	list size */
    if(rank == 0){
        if (0 > (list_size = read_input(input_name, &values))) return 2;
		extras = list_size%size;

		/* inserting zeroes to make the sort evenly divisable by 
		processors */
		if(extras!=0){
			list_size += extras;
			printf("---- OBS NOT DIVIDABLE BY %d --- extras=%d\n",
				size, extras);
			printf("placing %d zeores at end\n",extras);
			values = realloc(values, (list_size)*sizeof(int));
			for(int i=1;i<=extras;i++)	values[list_size-i]=0;
		}
		timing = (double*)malloc(size*sizeof(double));

		/* Pringting the input */
        printf("list size = %d\nfirst 32 elements in list: ",list_size);
        for(int i=0;i<32;i++){
            printf("%d ",values[i]);
        }   printf("\n");

		/* local sort test */
        qsort(values, list_size, sizeof(int), cmpfunc);

        /* printing the 32 first output values */
        printf("list size = %d\nfirst 32 elements in sorted list: ",list_size);
        for(int i=0;i<list_size;i++){
            printf("%d ",values[i]);
        }   printf("\n--------\n\n");
    }
	MPI_Bcast(&list_size, 1, MPI_INT, 0, hyperCube);
	MPI_Bcast(&extras, 1, MPI_INT, 0, hyperCube);
	local_size = list_size/size; // number of values per processor
	printf("..Broadcasting done!\n");

    //Defining the number of vectors to be sent (1 for every core) and the displacement which is different for each core
    int displacement[size];
	int count_send[size];
    for(int i =0; i<size; i++){
        displacement[i] = i*local_size;
		count_send[i] = local_size;
    }

    //Scattering A and B into local block matrices
	printf("Scattering local_values..\n");
	local_values = (int*)malloc(local_size*sizeof(int));
    MPI_Scatterv(values, count_send, displacement, MPI_INT, local_values, 
		local_size, MPI_INT, 0, hyperCube);
	printf("..Scattering done\n\n");

	if(rank==1){
		printf("local matrix in rank %d looks like: ", 1);
		for(int i=0;i<local_size;i++){
			printf("%d ", local_values[i]);
		}
		printf("\n-------\n\n");
	}
    


    /* writing the output */
    int output_status = 1; // 1 if active output. 0 inactive
	printf("new pointer skipp %d first values\n",extras);
	if(extras!=0){
		int *out = (&values)[extras];
		printf("after removeing first values\n");
	}
	
    if(rank == 0 && output_status==1){
		if (0 != write_output(output_name, values, list_size)) {
			return 2;
		}
		printf("after write output\n");
    }
    if(rank==0)
        free(values);
    MPI_Comm_free(&hyperCube);
    MPI_Finalize();
    return 0;
}




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

int write_output(const char *file_name, int *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%d ", output[i])) {
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