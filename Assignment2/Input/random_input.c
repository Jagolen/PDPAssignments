
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int write_output(const char *file_name, int *C, int matrix_size){
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < matrix_size * matrix_size * 2 + 1; i++) {
		if (0 > fprintf(file, "%d ", C[i])) {
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


int main(int argc, char **argv){
    
    if (3 != argc) {
        printf("Usage: matrix size N, output_file\n");
        return 1;
	}
    int N = atoi(argv[1]);
    char *output_name = argv[2];
    int matrix_size=N*N;
    int *matrix = (int*)malloc(matrix_size*2*sizeof(int));
    for(int i=1;i<matrix_size*2+1;i++){
        matrix[i]=rand()%100;
    }
    matrix[0]=N;
    write_output(output_name, matrix, N);

    return 1;
}