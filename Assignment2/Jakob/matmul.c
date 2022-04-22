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
    int period[2] = {1, 1};
    int size_per_dim[2], position[2];

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

    printf("My rank is %d and my coords are (%d, %d). \
Neighbors: North: %d, South: %d, East: %d, West: %d\n",
        rank, position[0], position[1],north, south,east,west);


    MPI_Finalize();

    return 0;
}