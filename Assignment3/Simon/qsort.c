#include "qsort.h"
int main(int argc, char **argv){
    
    if (3 != argc) {
        printf("Usage: input_file output_file\n");
        return 1;
	}
	char *input_name = argv[1];
	char *output_name = argv[2];
    double *timing;
    int rank, size, list_size;
    int *values;

    /* reading input where values is the list and list_size is the list size */
    if(rank == 0){
        if (0 > (list_size = read_input(input_name, &values))) return 2;
		timing = (double*)malloc(size*sizeof(double));
    }
    printf("list size = %d\nfirst 10 elements in list: ",list_size);
    for(int i=0;i<10;i++){
        printf("%d ",values[i]);
    }   printf("\n");


    /* writing the output */
    int output_status = 1; // 1 if active output. 0 inactive
    if(rank == 0 && output_status==1){
		if (0 != write_output(output_name, values, list_size)) {
			return 2;
		}
    }

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
	if (NULL == (*values = malloc(num_values * sizeof(double)))) {
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