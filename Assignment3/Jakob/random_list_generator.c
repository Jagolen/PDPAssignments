#include <stdio.h>
#include <stdlib.h>

int write_output(const char *filename, int *C, int leangth){
    FILE *file;
    if (NULL == (file = fopen(filename, "w"))){
        perror("could not open output file");
        return -1;
    }
    for (int i=0;i<=leangth;i++){
        if(0>fprintf(file, "%d ", C[i])){
            perror("Could not write output file");
        }
    }
    if (0>fprintf(file, "\n"))  perror("Could not write output file");
    if (0!=fclose(file)) perror("Could not close output file");
    return 0;
}

int main(int argc, char **argv){
    if(3!=argc){
        printf("Usage: list leangth N, output file\n");
        return 1;
    }
    int N = atoi(argv[1]);
    char *output_name = argv[2];
    if(N<=1){
        printf("N must be large than 1\n");
        return 1;
    }
    int *vector = (int *)malloc((N+1)*sizeof(int));
    for(int i=1;i<=N;i++) vector[i]=rand()%100 - 50;
    vector[0]=N;
    write_output(output_name, vector, N);
    return 0;
}