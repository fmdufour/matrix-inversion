#include "util.h"
#include "matrix.h"

int main(int argc, char* argv[]){        
    configuration config;        
    double* matrix;
    double* invMatrix;
    double* Y;
    double *B;
    double* idMatrix;    
    pivotsRecord* pRecord;    
    srand(20172);
    LU* lu;    
    
    config = readConfiguration(argc, argv);    
        
    if(config.inputFile == NULL)
        matrix = generateSquareRandomMatrix(config.randomMatrixSize);
    else
        matrix = readMatrixFromFile(config);

    pRecord = malloc(sizeof(pivotsRecord));
    pRecord->pivots = malloc(sizeof(pivot)* N); 
    pRecord->count = 0;

    lu = luDecomposition(matrix, pRecord, N);
                
    idMatrix = getIdentityMatrix(N);    
    pivotMatrix(idMatrix, pRecord);

    Y = forwardSubstitution(lu->L, idMatrix, N, N);        
    
    invMatrix = backwardSubstitution(lu->U, Y, N, N);       

    B = multiplyMatrix(matrix, invMatrix, N);
    
    printf("-----A------\n");
    printMatrix(matrix, N);    
    printf("-----L------\n");
    printMatrix(lu->L, N);
    printf("-----U------\n");
    printMatrix(lu->U, N);    
    printf("-----Y------\n");    
    printMatrix(Y, N);
    printf("---A-1----\n");
    printMatrix(invMatrix, N);
    printf("-----AxA-1------\n");    
    printMatrix(B, N);    

    free(Y);
    free(idMatrix);
    free(matrix);
    free(lu->L);
    free(lu->U);
    free(lu);

    return 0;
}