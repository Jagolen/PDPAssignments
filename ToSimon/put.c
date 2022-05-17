#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
    int rank,size;
    double *output, *timings;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Win win;
    MPI_Alloc_mem(size*1000*sizeof(double), MPI_INFO_NULL, &output);
    MPI_Win_create(output, size*1000*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    printf("Starting timing\n");
    timings = (double*)malloc(1000*sizeof(double));
    
    for(int i=0; i<1000; i++){
        double timing = MPI_Wtime();
        long int a = 1;
        for(long int j=0; j<1000000; j++){
            a++;
        }
        timing = timing - MPI_Wtime();
        timings[i] = timing;
    }
    printf("Timing done!\n");

    for(int i = 0; i<1000; i++){
        MPI_Win_fence(MPI_MODE_NOPRECEDE,win);
        MPI_Put(timings,1000,MPI_DOUBLE,0,1000*rank,1000,MPI_DOUBLE, win);
        MPI_Win_fence(MPI_MODE_NOSUCCEED,win);
    }
    
    MPI_Finalize();
    return 0;
}