#     include <mpi.h>
#     include <stdio.h>

/* Run with FOUR processes */
void main(int argc, char *argv[]) 
{
        int rank, size;
        MPI_Status status;
        struct {
              int x;  int y;   int z;
               } point;
        MPI_Datatype ptype;
        MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
        MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
	if(size!=4){
           printf("Run with four processes!\n");             /* Print a message              */
           MPI_Finalize();
	   return; 
	}
        MPI_Type_contiguous(3,MPI_INT,&ptype);
        MPI_Type_commit(&ptype);
	point.x=0+rank; point.y=1+rank; point.z=2+rank;
        printf("P:%d contents of x,y,z is (%d,%d,%d) \n",rank,point.x,point.y,point.z);
	
        if(rank==3){
           point.x=15; point.y=23; point.z=6;
           MPI_Send(&point,1,ptype,1,52,MPI_COMM_WORLD);
        } else if(rank==1) {
        MPI_Recv(&point,1,ptype,3,52,MPI_COMM_WORLD,&status);
        printf("P:%d received coords are (%d,%d,%d) \n",rank,point.x,point.y,point.z);
     }
     MPI_Finalize();
}
     
