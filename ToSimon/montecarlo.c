/* 
Monte Carlo with the Stochastic Simulation Algorithm
By: Jakob Gölén
*/

#include "montecarlo.h"

int main(int argc, char **argv){
    if (3 != argc) {
        printf("Usage: num_experiments output_file\n");
        return 1;
	}

	//Parameters

	//Inputs
    int num_experiments = atoi(argv[1]);
    char *output_name = argv[2];

	//To store timings
	double local_time[5] = {0, 0, 0, 0, 0};
	double *times;
	double start_final, start_sub, end_final;

	//To determine bin edges
	int local_min, local_max, global_min, global_max;
	int bins[21];
	int bincounts[20] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	//Given constant
	const double T = 100;

	//Parameters used in the time loop
	double u1, u2, a0, tau, tt, r_data;
	int r, quarter_done, half_done, three_quarter_done, all_done;

	//State and result vectors
	double *w;
	int *x, *p;
	int *result;
	int *local_result;

	//For MPI
	int rank, size;

	//Initializing MPI
	MPI_Init(&argc, &argv);

	//Getting rank and size
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//number of experiments per process
	int n = num_experiments/size;
	if(num_experiments%size != 0){
		if(rank == 0) printf("Number of experiments (%d) is not divisible by number of processes (%d)\n", num_experiments, size);
		MPI_Abort(MPI_COMM_WORLD, -1);
		exit(-1);
	}

	//Initializing windows to put the result and timings to rank 0, and allocates the memory
	MPI_Win reswin;
	MPI_Win timewin;

	if(rank == 0){
		MPI_Alloc_mem(20*size*sizeof(int), MPI_INFO_NULL, &result);
		MPI_Win_create(result, 20*size*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &reswin);
		MPI_Alloc_mem(5*size*sizeof(double), MPI_INFO_NULL, &times);
		MPI_Win_create(times, 5*size*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &timewin);
	}
	else{
		MPI_Win_create(NULL, 0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &reswin);
		MPI_Win_create(NULL, 0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &timewin);
	}


	//Initializing x, w and p
	x = (int*)malloc(7*sizeof(int));
	w = (double*)malloc(15*sizeof(double));
	p = (int*)malloc(15*7*sizeof(int));

	//Allocating local result vector
	local_result = (int*)malloc(n*sizeof(int));

	//the initial values of x was given, the values are initialized
	const int x0[7] = {900, 900, 30, 330, 50, 270, 20};

	//P is a constant and most values are zero. The zero and non zero elements are initialized
	for(int i = 0; i<(15*7);i++) p[i] = 0;

	p[0] = 1; p[7] = -1; p[14] = -1; p[16] = 1; p[22] = 1; p[29] = -1; p[36] = -1; p[38] = 1; p[44] = -1; p[51] = -1; p[53] = 1; 
	p[59] = -1; p[66] = -1; p[68] = 1; p[74] = -1; p[81] = -1; p[83] = 1; p[89] = -1; p[91] = 1; p[97] = -1; p[104] = -1;

	//To generate different numbers every time the program is ran and gives different random numbers on the different processes
	srand(time(0)+rank);

	//Fences so the data has access to the windows
	MPI_Win_fence(MPI_MODE_NOPRECEDE, reswin);
	MPI_Win_fence(MPI_MODE_NOPRECEDE, timewin);

	start_final = MPI_Wtime();
	//Main loop
	for(int iter = 0; iter<n; iter++){
		//Time is reset and x is set to initial value
		tt = 0;
		for(int i = 0; i<7; i++) x[i] = x0[i];
		quarter_done = 0;
		half_done = 0;
		three_quarter_done = 0;
		all_done = 0;

		//Time loop for every iteration
		start_sub = MPI_Wtime();
		while(tt<T){

			//Calculating w using the prop function and then calculating the sum of w
			a0 = 0;
			prop(x,w);
			for(int i = 0; i<15;i++) a0 += w[i];

			/*
			Draw two random uniform numbers between 0 and 1. 0 is excluded since log 0 is undefined

			*/
			u1 = ((double)(rand()%(RAND_MAX-1)+1))/(RAND_MAX);
			u2 = ((double)(rand()%(RAND_MAX-1)+1))/(RAND_MAX);

			//Finding tau and r
			tau = -log(u1)/a0;
			r = 0;
			r_data = w[0];
			while(r_data < (a0*u2)){
				r++;
				r_data+=w[r];
			}

			//Updating x based on r
			for(int i = 0; i<7; i++) x[i] += p[7*r+i];

			//Updating the time
			tt += tau;

			//If a subinterval is passed, add to the mean time of that subinterval
			if(quarter_done == 0 && tt>25){
				local_time[0] += ((MPI_Wtime() - start_sub)/n);
				quarter_done = 1;
			}

			if(tt>50 && half_done == 0){
				local_time[1] += ((MPI_Wtime() - start_sub)/n);
				half_done = 1;
			}

			if(tt>75 && three_quarter_done == 0){
				local_time[2] += ((MPI_Wtime() - start_sub)/n);
				three_quarter_done = 1;
			}

		}
		local_time[3] += ((MPI_Wtime() - start_sub)/n);

		//Saving the x[0] element in the result vector
		local_result[iter] = x[0];

		//Finding min and max value of x[0]
		if(iter == 0){
			local_min = x[0];
			local_max = x[0];
		}

		if(x[0] < local_min) local_min = x[0];
		if(x[0] > local_max) local_max = x[0];

	}

	//The result has been collected, timer stops
	local_time[4] = MPI_Wtime() - start_final;

	//Sending the times to rank 0 with put
	MPI_Put(local_time, 5, MPI_DOUBLE, 0, rank*5, 5, MPI_DOUBLE, timewin);

	//Next, the data is prepared to be shown in a histogram

	//Getting the global min and max values
	MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);


	//All processes now have tha global max and min, determining the bin vector
	int bin_step_size = (global_max-global_min)/20;
	int current_bin = global_min;
	bins[0] = global_min; bins[20] = global_max;
	for(int i = 1; i<20; i++){
		current_bin+=bin_step_size;
		bins[i] = current_bin;
	}

	//Now the results are divided into the bins
	int chosen_bin;
	for(int i = 0; i<n; i++){
		chosen_bin = 0;
		while(local_result[i] > bins[chosen_bin+1]) chosen_bin++;
		bincounts[chosen_bin]++;
	}

	//Sending the bin counts to process 0
	MPI_Put(bincounts, 20, MPI_INT, 0, rank*20, 20, MPI_INT, reswin);


	//Putting fences so the data can be accessed
	MPI_Win_fence(MPI_MODE_NOSUCCEED, reswin);
	MPI_Win_fence(MPI_MODE_NOSUCCEED, timewin);

	//printing the times
	if(rank == 0){
		printf("Average times: \n");
		printf("Rank\t T=25 \t\t T=50 \t\t T=75 \t\t T=100 \n");
		for(int i = 0; i<size; i++) printf("%d \t %lf \t %lf \t %lf \t %lf \n", i, times[5*i], times[5*i+1], times[5*i+2], times[5*i+3]);

		//Print the highest execution time
		double max_time = times[4];
		for(int i = 1; i<size; i++){
			if(max_time < times[5*i+4]) max_time = times[5*i+4];
		}
		printf("Highest time: %lf\n", max_time);

		//Merging the bins
		int final_bins[20] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		for(int i = 0; i<20*size; i++){
			final_bins[i%20] += result[i];
		}



		//Write the output file
        if (0 != write_output(output_name, result, 20, bins)) {
            return 2;
		}
	}

	//Freeing memory
	free(x);
	free(w);
	free(p);


	//Freeing the windows
	MPI_Win_free(&reswin);
	MPI_Win_free(&timewin);
	if(rank == 0){
		MPI_Free_mem(result);
		MPI_Free_mem(times);
	}

	//Finalizing MPI and end the program
	MPI_Finalize();
    return 0;
}

//Prop function used to find w
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


//Function for writing the output
int write_output(const char *file_name, int *result, int num_bins, int *bins) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_bins; i++) {
		if (0 > fprintf(file, "%d ", result[i])) {
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "\n")) {
		perror("Couldn't write to output file");
	}
	
	for(int i = 0; i<20; i++){
		if (0 > fprintf(file, "%d ", bins[i])) {
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