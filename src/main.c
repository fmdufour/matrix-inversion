#include "util.h"
#include "matrix.h"

int main(int argc, char* argv[]){        
    configuration config;        
    double* matrix;
    double* identityColumn;
    double* invMatrix;
    double* Y;
    int i, j;
    srand(20172);
    LU* lu;
    
    config = readConfiguration(argc, argv);    

    printf("Configurações recebidas:\n");    
    printf("Arquivo de Entrada: %s\n", config.inputFilePath);
    printf("Arquivo de Saida: %s\n", config.outputFilePath);
    printf("Tam Matriz aleatoria: %d\n", config.randomMatrixSize);
    printf("Numero de Iteracoes Refinamento: %d\n", config.iterationCount);
        
    if(config.inputFile == NULL)
        matrix = generateSquareRandomMatrix(config.randomMatrixSize);
    else
        matrix = readMatrixFromFile(config);
    
    //printMatrix(matrix, N);

    lu = luDecomposition(matrix, N);
    
    Y = allocateMatrix(N);
    invMatrix = allocateMatrix(N);

    for(i = 0; i < N; i ++)
        getVal(lu->L, i, i) = 1.0;

    printf("-----L------\n");
    printMatrix(lu->L, N);
    printf("-----U------\n");
    printMatrix(lu->U, N);    

    for(i = 0; i < N; i ++){                
        forwardSubstitution(lu->L, Y, i, N);
        //for(j = 0; j < N; j++)
            //printf("%f ", identityColumn[j]);
    }
    printf("-----Y------\n");    
    printMatrix(Y, N);
    for(i = 0; i < N; i ++){                
        backwardSubstitution(lu->U, invMatrix, Y, i, N);        
        //for(j = 0; j < N; j++)
            //printf("%f ", identityColumn[j]);
    }      

    // printf("--L--\n");
    // printMatrix(lu->L, N);
    // printf("--L--\n");    
    printf("---A-1----\n");
    printMatrix(invMatrix, N);
    // printf("--A--\n");    
    // printMatrix(matrix, N);
    // printf("--L--\n");
    // printMatrix(lu->L, N);
    // printf("--U--\n");
    // printMatrix(lu->U, N);

    return 0;
}