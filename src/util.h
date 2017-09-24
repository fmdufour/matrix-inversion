#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

int N;
#define getVal(m, i, j) (m[j + i * N])

typedef struct configuration{
    char* inputFilePath;
    FILE* inputFile;
    FILE* outputFile;    
    char* outputFilePath;    
    unsigned int randomMatrixSize;
    int iterationCount;
} configuration;

double* getIdentityColumn(unsigned int iCol, unsigned int n);
double* allocateMatrix(unsigned int n);
double *generateSquareRandomMatrix(unsigned int n);
double *readMatrixFromFile(configuration config);
double timestamp(void);
void printMatrix(double* matrix, unsigned int n);
void printErrorExit(char* msg);
configuration readConfiguration(int arg, char *argv[]);