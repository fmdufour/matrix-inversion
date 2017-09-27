#include "util.h"
#include "matrix.h"

/**
 * Autores
 * 
 * @autor Fernando Medeiros Dufour GRR20140513
 *   
 * @autor Carolina Aparecida de Lara Moraes GRR20111353           
 */


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
    double startTime, startTime2;
    double luDecompTime = 0.0;    
    double totalResidueTime = 0.0;
    double totalLsSolveTime = 0.0;
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

    startTime = timestamp();

    lu = luDecomposition(matrix, pRecord, N);

    luDecompTime = diffInSeconds(startTime);

    pivotMatrix(matrix, pRecord);
                
    idMatrixPivoted = getIdentityMatrix(N);

    Y = forwardSubstitution(lu->L, idMatrixPivoted, N, N);        
    
    invMatrix = backwardSubstitution(lu->U, Y, N, N);       
    
    idMatrix = getIdentityMatrix(N);            

    if(config.outputFile == NULL)
        printf("#\n");
    else
        fprintf(config.outputFile, "#\n");    

    for(i = 0; i < config.iterationCount; i++){           
        startTime = timestamp();

        B = multiplyMatrix(matrix, invMatrix, N);

        startTime2 = timestamp();

        R = getResidue(B, idMatrix, &residue, N);                        

        totalResidueTime += diffInSeconds(startTime2);
        
        Y = forwardSubstitution(lu->L, R, N, N);

        Delta = backwardSubstitution(lu->U, Y, N, N);

        AddMatrix(invMatrix, Delta, N);

        totalLsSolveTime += diffInSeconds(startTime);

        if(config.outputFile == NULL)
            printf("# iter %d: %.17f\n", (i+1), residue);
        else
            fprintf(config.outputFile, "# iter %d: %.17f\n", (i+1), residue);

    }        
        
    pivotMatrix(invMatrix, pRecord);    

    if(config.outputFile == NULL){
        printf("# Tempo LU: %.17f\n", luDecompTime);
        printf("# Tempo iter: %.17f\n", totalLsSolveTime/config.iterationCount);
        printf("# Tempo residuo: %.17f\n", totalResidueTime/config.iterationCount);
        printf("#\n");
        printMatrix(invMatrix, N);        
    }
    else{
        fprintf(config.outputFile, "# Tempo LU: %.17f\n", luDecompTime);
        fprintf(config.outputFile, "# Tempo iter: %.17f\n", totalLsSolveTime/config.iterationCount);
        fprintf(config.outputFile, "# Tempo residuo: %.17f\n", totalResidueTime/config.iterationCount);
        fprintf(config.outputFile, "#\n");
        printMatrixToFile(invMatrix, N, config.outputFile);
    }

    return 0;
}
