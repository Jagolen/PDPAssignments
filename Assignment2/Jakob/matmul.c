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
    int north, south, east, west;
    const int *period = {1, 1};
    const int *size_per_dim, *position;

    MPI_Init(&argc, &argv);
    
    
    //Creating a cartesian topology and get the ranks
    MPI_Comm cart;


    MPI_Finalize();

    return 0;
}