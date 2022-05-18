#include "mcssa.h"



int main(int argc, char **argv){
    if(argc != 3){
        printf("Usage: N experiments, T time\n");
        return 1;
    }
    int x[7] = {0,0,0,0,0,0,0};
    double times;
    double *w;
    int rank, size, n;
    const int N = atoi(argv[1]);
    const double T = atoi(argv[2]);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD , &size);
    if(N%size!=0){
        printf("N (%d) not evenly divisible by size (%d). exiting program\n",
            N, size);
        return 1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD , &rank);
    
    if(rank==0){
        int n = N/size;
        printf("local experiments = %d\n",n);
    }
    MPI_Bcast(&n ,1 ,MPI_INT ,0 , MPI_COMM_WORLD);
    double t=0, a0, u1, u2, dt, less, more;
    int i, r;
    srand(rank*100);
 
    w = (double *)malloc(15*sizeof(double));
    while(t<T){
        prop(x, w);
        a0=0;
        for(i=0;i<15;i++) a0 += w[i];
        //printf("a0=%f\n",a0);
        u1 = (double)rand()/RAND_MAX;
        u2 = (double)rand()/RAND_MAX;
        dt = -log(u1)/a0;
        less = more = 0;
        /*while(!(less<a0*u2 && more>a0*u2)){
            if(r>0) less += w[r-1];
            more += w[r]
            r++
        }*/
        if(rank==0){
        int *P_vals = (int*)malloc(21*sizeof(int));
        int *colms = (int*)malloc(7*sizeof(int));
        int *rows = (int*)malloc(21*sizeof(int));
        int *P = (int*)malloc(7*sizeof(int));
        P_matrix(P_vals, colms, rows);
            
        }

        t++;
    }
    



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

void P_matrix(int *P_vals, int *colms, int *rows){
    /* CSC sparse matrix */
    int rwws[21] = {0, 1, 2, 13, 3, 4, 5, 2, 6, 7, 5, 8, 9
                    , 7, 10, 11, 9, 12, 11, 13, 14};
    int valuess = {1, -1, -1, 1, 1, -1, -1, 1, -1, -1, 1, 
                    -1, -1, 1, -1, -1, 1, -1, 1, -1, -1};
    for(int i=0;i<21;i++){
        /* setting P values to allocated memmory */
        if(i%3==0) P_vals[i] = values[i];
        printf("%d ",P_vals[i]);
        /* setting row values to allocated memmory */
        rows[i] = rwws[i];
    } printf("\n");
    
    for(int j=1;j<6;j++){
        colms[j] = 1+j*3;
        printf("%d ",colms[j]);
    } printf("\n");;
    colms[0]=0; colms[6]=18; colms[7]=21;
    
}

void row_from_sparse_converter(int r, int *P_vals, 
    int *colms, *int rows, int *P){
    for(int k= 0;k<7;k++){
        P[k]=0;
    }
    int j;
    for(int i=0;i<21;i++){
        if(rows[i]==r){
            j=0;
            while(colms[j]<i+1){
                j++;
                if(colms[j]>i+1) P[j]=P_vals[i];
            }
        }
    }
}
