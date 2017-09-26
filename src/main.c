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
    double time;
    pivotsRecord* pRecord;    
    srand(20172);
    LU* lu;    
    
    config = readConfiguration(argc, argv);    
        
    if(config.randomMatrixSize > 0 && config.inputFile == NULL)
        matrix = generateSquareRandomMatrix(config.randomMatrixSize);
    else if(config.inputFile != NULL)
        matrix = readMatrixFromFile(config);
    else
        matrix = readMatrixFromStdIn();

    pRecord = malloc(sizeof(pivotsRecord));
    pRecord->pivots = malloc(sizeof(pivot)* N); 
    pRecord->count = 0;

    lu = luDecomposition(matrix, pRecord, N);

    pivotMatrix(matrix, pRecord);
                
    idMatrixPivoted = getIdentityMatrix(N);

    Y = forwardSubstitution(lu->L, idMatrixPivoted, N, N);        
    
    invMatrix = backwardSubstitution(lu->U, Y, N, N);       
    
    idMatrix = getIdentityMatrix(N);            

    // printMatrix(lu->L, N);
    // printf("\n");
    // printMatrix(lu->U, N);
    // printf("#\n");
    //double startTime = timestamp();
    //Mprintf("%.17f\n", teste);    
//    printMatrix(invMatrix, N);

    for(i = 0; i < config.iterationCount; i++){
        //pivotMatrix(invMatrix, pRecord);

        B = multiplyMatrix(matrix, invMatrix, N);

        R = getResidue(B, idMatrix, &residue, N);                        

        
        Y = forwardSubstitution(lu->L, R, N, N);

        Delta = backwardSubstitution(lu->U, Y, N, N);

        AddMatrix(invMatrix, Delta, N);

        printf("# iter %d: %.17f\n", (i+1), residue);

    }    
        
    // printMatrix(B, N);
    // printf("%d\n", N);    
    //printMatrix(invMatrix, N);    

    //printf("%.17f", fabs(teste - timestamp()/1000));

    return 0;
}