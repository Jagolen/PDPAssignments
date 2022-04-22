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
    int *period = {1, 1};
    int *size_per_dim[2], *position[2];

    //Initializing MPI and get the number of processing elements
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    
    /*Creating a cartesian topology and get the ranks a 2d square topology is assumed
    so processes per dimension is just sqrt(size)*/
    MPI_Comm cart;
    size_per_dim[0] = sqrt(size);
    size_per_dim[1] = sqrt(size);


    MPI_Finalize();

    return 0;
}