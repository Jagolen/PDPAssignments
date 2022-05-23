#include "mcssa.h"



int main(int argc, char **argv){
    if(!(argc==2 || argc==3)){
        printf("Usage: N experiments and optionally output name.\n");
        return 1;
    }
    const int xx[7] = {900,900,30,330,50,270,20};
    const int N = atoi(argv[1]);
    int rank, size, n, i, r;
    int *P_vals, *results, *x;
    int output_status = 0;
    const char *output_name;
    if(argc == 3){output_name=argv[2]; output_status=1;}
    const double T = 100;
    double a0, u1, u2, dt, less, more, t, start_time, 
        start_rank_time, rank_time;
    double *w, *rank_times;
    double times[8];


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD , &size);
    if(N%size!=0){
        printf("N (%d) not evenly divisible by size (%d). exiting program\n",
            N, size+10);
        return 1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD , &rank);
    MPI_Win win;
    MPI_Win time_win;

    /* setting local size and broadcasting it to all processors */
    if(rank==0){
        n = N/size;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    srand(1003*rank);
    

    /* allocating and creating windows */
    if(rank==0){
        /* Window for gathering histogram result data */
        MPI_Alloc_mem(N*sizeof(int), MPI_INFO_NULL, &results);
        MPI_Win_create(results , N*sizeof(int), sizeof(int),
            MPI_INFO_NULL , MPI_COMM_WORLD , &win);
        
        /* window for gathering total rank times */
        MPI_Alloc_mem(size*sizeof(double), MPI_INFO_NULL, &rank_times);
        MPI_Win_create(rank_times, size*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD, &time_win);
    }
    else{
        MPI_Win_create(NULL, 0, sizeof(int),
            MPI_INFO_NULL, MPI_COMM_WORLD, &win);
        MPI_Win_create(NULL, 0, sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD, &time_win);
    }
    P_vals = (int*)malloc(15*7*sizeof(int));   //allocating P-matrix
    P_matrix(P_vals); // Setting values of P matrix
    w = (double *)malloc(15*sizeof(double));
    x = (int*)malloc(7*sizeof(int));

    /* fences to allow the windows to be accessed */
    MPI_Win_fence(MPI_MODE_NOPRECEDE, win);
    MPI_Win_fence(MPI_MODE_NOPRECEDE, time_win);

    for(i=0;i<4;i++) times[i]=0;  // timing values to 0;


    start_rank_time = MPI_Wtime();
    /* Main loop n local iterations over time T with random time steps */
    for(int iter=0;iter<n;iter++){
        t=0;
        for(i=0;i<7;i++) x[i]=xx[i]; // x to x0
        for(i=4;i<8;i++) times[i]=0;    // time has not been uppdated
        start_time = MPI_Wtime();
        while(t<T){
            prop(x, w);
            a0=0;
            for(i=0;i<15;i++) a0 += w[i];
            u1 = ((double)(rand()%1000000)+1)/(1000001);
            u2 = ((double)(rand()%1000000)+1)/(1000001);
            dt = -log(u1)/a0;
            more = 0;
            less = 0;
            r = 0;
            
            /*  finding w(i) sums between a0*u2 */

            while(!(less<a0*u2 && more>=a0*u2)){
                more+=w[r];
                if(r>0) less+=w[r-1];
                if(r>14){
                    printf("%f %f %f\n\n",less, a0*u2, more);
                    printf("stuck in local while\n");
                    return 1;
                }
                r++;
            }
            /* update x vector */
            for(i=0;i<7;i++){
                x[i] += P_vals[(r-1)*7+i];
            }
            t+=dt; //Update time
            /* getting mean time after time variable 25, 50, 75 and 100. */
            if(times[4]==0 && t>25){    //time is over 25 and not been updated
                times[0]+=(MPI_Wtime()-start_time)/n;   //time is set
                times[4]=1; // time for 25 has been updated
            }
            else if(t>50 && times[5]==0){
                times[1]+=(MPI_Wtime()-start_time)/n;
                times[5]=1;
            }
            else if(t>75 && times[6]==0){
                times[2]+=(MPI_Wtime()-start_time)/n;
                times[6]=1;
            }
            else if(t>75 && times[7]==0){
                times[3]+=(MPI_Wtime()-start_time)/n;
                times[7]=1;
            }
        }
        MPI_Put(&x[0], 1, MPI_INT, 0, rank*n+iter, 1, 
                MPI_INT, win);   
    }
    /* gathering total rank times to rank 0 */
    rank_time = MPI_Wtime()-start_rank_time;
    MPI_Put(&rank_time, 1, MPI_DOUBLE, 0, rank, 1, MPI_DOUBLE,
             time_win);

    /* makes sure everything is sent to rank 0*/
    MPI_Win_fence(MPI_MODE_NOSUCCEED,win);
    MPI_Win_fence(MPI_MODE_NOSUCCEED,time_win);
    MPI_Barrier(MPI_COMM_WORLD);

    /*  Printing out the times in the order of 25, 50, 75 and 100 
        along the rows and the ranks along the column axis */
    printf("rank %d: %.3f %.3f %.3f %.3f (ms)\n",rank, times[0]*1000, times[1]*1000
        , times[2]*1000, times[3]*1000);
    /* printing total time of ranks */    
    if(rank==0) {
        double max_time =0;
        for(int i=0;i<size;i++){
            if(rank_times[i]>max_time) max_time=rank_times[i];
        }
        printf("slowest total rank time = %.3fs\n",max_time);
    }

    

        /* writing the output file */
    if(output_status==1 && rank == 0){
		if (0 != write_output(output_name, results, N)) {
            printf("could not write output\n");
			return 2;
		}
    }
    /* freeing memmory and closing windows */
    if(rank==0) {MPI_Free_mem(results); MPI_Free_mem(rank_times);}
    MPI_Win_free(&win); MPI_Win_free(&time_win)  
    free(w); free(P_vals), free(x);
    MPI_Finalize();
    return 0;
}



/* properties to set the values of w dependent on x */
void prop(int *x, double *w) {
	// Birth number, humans
	const double LAMBDA_H = 20;
	// Birth number, mosquitoes
	const double LAMBDA_M = 0.5;
	// Biting rate of mosquitoes
	const double B = 0.075;
	/* Probability that a bite by an infectious mosquito results in transmission
	   of disease to human*/
	const double BETA_H = 0.3;
	/* Probability that a bite results in transmission of parasite to a
	   susceptible mosquito*/
	const double BETA_M = 0.5;
	// Human mortality rate
	const double MU_H = 0.015;
	// Mosquito mortality rate
	const double MU_M = 0.02;
	// Disease induced death rate, humans
	const double DELTA_H = 0.05;
	// Disease induced death rate, mosquitoes
	const double DELTA_M = 0.15;
	// Rate of progression from exposed to infectious state, humans
	const double ALFA_H = 0.6;
	// Rate of progression from exposed to infectious state, mosquitoes
	const double ALFA_M = 0.6;
	// Recovery rate, humans
	const double R = 0.05;
	// Loss of immunity rate, humans
	const double OMEGA = 0.02;
	/* Proportion of an antibody produced by human in response to the incidence
	   of infection caused by mosquito. */
	const double NU_H = 0.5;
	/* Proportion of an antibody produced by mosquito in response to the
	   incidence of infection caused by human. */
	const double NU_M = 0.15;

	w[0] = LAMBDA_H;
	w[1] = MU_H * x[0];
	w[2] = (B * BETA_H * x[0] * x[5]) / (1 + NU_H * x[5]);
	w[3] = LAMBDA_M;
	w[4] = MU_M * x[1];
	w[5] = (B * BETA_M * x[1]*x[4]) / (1 + NU_M * x[4]);
	w[6] = MU_H * x[2];
	w[7] = ALFA_H * x[2];
	w[8] = MU_M * x[3];
	w[9] = ALFA_M * x[3];
	w[10] = (MU_H + DELTA_H) * x[4];
	w[11] = R * x[4];
	w[12] = (MU_M + DELTA_M) * x[5];
	w[13] = OMEGA * x[6];
	w[14] = MU_H * x[6];
}

void P_matrix(int *P_vals){
    /* CSC sparse matrix */
    int num_cols=7;
    for(int i=0;i<7*15;i++) P_vals[i]=0;
    P_vals[0]=1; P_vals[7*1]=-1; P_vals[2*7]=-1;
    P_vals[2*7+2]=1; P_vals[3*7+1]=1; P_vals[4*7+1] = -1;
    P_vals[5*7+1]=-1; P_vals[5*7+3]; P_vals[6*7+2]=-1;
    P_vals[7*7+2]=-1; P_vals[7*7+4]=1; P_vals[8*7+3]=-1;
    P_vals[9*7+3]=-1; P_vals[9*7+5]=1; P_vals[10*7+4]=-1;
    P_vals[11*7+4]=-1; P_vals[11*7+6]=1; P_vals[12*7+5]=-1;
    P_vals[13*7]=1; P_vals[13*7+6]=-1; P_vals[14*7+6]=-1;
}

int write_output(const char *file_name, int *results, int N) {
	FILE *file;
    /* Create a file with file_name and ready to write */
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
    /* write values from results list to file */
	for (int i = 0; i < N; i++) {
		if (0 > fprintf(file, "%d ", results[i])) {
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