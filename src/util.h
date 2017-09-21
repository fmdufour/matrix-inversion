#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

typedef struct configuration{
    char* inputFile;
    char* outputFile;
    unsigned int randomMatrixSize;
    int iterationCount;
} configuration;

double** generateSquareRandomMatrix( unsigned int n );
double timestamp(void);
void printMatrix(double** matrix, unsigned int n);
void printErrorExit(char* msg);