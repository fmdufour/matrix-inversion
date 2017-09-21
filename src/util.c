#include "util.h"

double** generateSquareRandomMatrix( unsigned int n )
{
    int i, j;
    double** mat = NULL;
    /* return NULL if memory allocation fails */
    if (!(mat = (double **) malloc(n*sizeof(double*))))
        return (NULL);
    for(i = 0; i < n; i ++) 
        if(!(mat[i] = (double*) malloc(n * sizeof(double*))))
            return NULL;

    /* generate a randomly initialized matrix in row-major order */    
    double invRandMax = 1.0/(double)RAND_MAX;

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++)
            mat[i][j] = (double)rand() * invRandMax;
    }
    
    return mat;
}



double timestamp(void){
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

void printMatrix(double** matrix, unsigned int n){
    int i, j;

    for(i = 0; i < n; i ++){
        for(j = 0; j < n; j++){
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }            
}

void printErrorExit(char* msg){
    fprintf(stderr, msg);
    exit(1);
}