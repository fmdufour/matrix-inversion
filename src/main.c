#include "util.h"
#include "matrix.h"


int main(int argc, char* argv[]){        
    configuration config;  
    int i = 0;      
    double* matrix;
    double* invMatrix;
    double* Y;
    double* B;
    double* R;
    double* Delta;
    double* idMatrixPivoted;    
    double* idMatrix;    
    double residue;
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
                
    idMatrixPivoted = getIdentityMatrix(N);
    pivotMatrix(idMatrixPivoted, pRecord);

    Y = forwardSubstitution(lu->L, idMatrixPivoted, N, N);        
    
    invMatrix = backwardSubstitution(lu->U, Y, N, N);       

    B = multiplyMatrix(matrix, invMatrix, N);
    idMatrix = getIdentityMatrix(N);        
    
    printf("#\n");

    for(i = 0; i < config.iterationCount; i++){

        R = getResidue(B, idMatrix, &residue, N);                        

        pivotMatrix(R, pRecord);
        
        Y = forwardSubstitution(lu->L, R, N, N);

        Delta = backwardSubstitution(lu->U, Y, N, N);

        AddMatrix(invMatrix, Delta, N);

        free(B);
        
        B = multiplyMatrix(matrix, invMatrix, N);
        
        printf("# iter %d: %.17f\n", (i+1), residue);

    }    
        
    printf("%d\n", N);    
    printMatrix(invMatrix, N);    

    free(Y);
    free(idMatrix);
    free(matrix);
    free(lu->L);
    free(lu->U);
    free(lu);

    return 0;
}