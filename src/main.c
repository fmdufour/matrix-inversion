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
    double* idMatrix;    
    double residue;
    double startTime, startTime2;
    double luDecompTime = 0.0;    
    double totalResidueTime = 0.0;
    double totalLsSolveTime = 0.0;
    pivotsRecord* pRecord;    
    srand(20172);
    LU* lu;    
    
    //Le do stdin as configuracoes e preenche a estrutura de dados
    config = readConfiguration(argc, argv);    
    
    //recupera a matriz em forma de vetor conforme a entrada por linha de comando
    if(config.randomMatrixSize > 0 && config.inputFile == NULL)
        matrix = generateSquareRandomMatrix(config.randomMatrixSize);
    else if(config.inputFile != NULL)
        matrix = readMatrixFromFile(config);
    else
        matrix = readMatrixFromStdIn();

    //aloca a estrutura de dados a qual guarda o historico de pivots
    pRecord = malloc(sizeof(pivotsRecord));
    pRecord->pivots = malloc(sizeof(pivot)* N); 
    pRecord->count = 0;

    startTime = timestamp();

    //fatora a matriz em L e U
    lu = luDecomposition(matrix, pRecord, N);    

    luDecompTime = diffInSeconds(startTime);

    //pivotei a matriz original para que a multiplicacao pela inversa gere a identidade corretamente
    pivotMatrix(matrix, pRecord);    
                
    //gera uma matriz identidade    
    idMatrix = getIdentityMatrix(N);    

    //soluciona a equacao da forma Ly=B
    Y = forwardSubstitution(lu->L, idMatrix, N, N);        
    //soluciona a equacao da forma Ux=Y
    invMatrix = backwardSubstitution(lu->U, Y, N, N);       

    if(config.outputFile == NULL)
        printf("#\n");
    else
        fprintf(config.outputFile, "#\n");    

    for(i = 0; i < config.iterationCount; i++){           
        startTime = timestamp();
        //realiza a multiplicacao da matriz por sua inversa
        B = multiplyMatrix(matrix, invMatrix, N);

        startTime2 = timestamp();
        // gera a matriz com os erros para cada x e tambem calcula residuo
        R = getResidue(B, idMatrix, &residue, N);                        

        totalResidueTime += diffInSeconds(startTime2);
        //resolve o sistema para o residuo
        Y = forwardSubstitution(lu->L, R, N, N);
        //resolve o sistema para gera o delta
        Delta = backwardSubstitution(lu->U, Y, N, N);
        // adiciona o delta na matriz inversa
        AddMatrix(invMatrix, Delta, N);

        totalLsSolveTime += diffInSeconds(startTime);

        if(config.outputFile == NULL)
            printf("# iter %d: %.17g\n", (i+1), residue);
        else
            fprintf(config.outputFile, "# iter %d: %.17g\n", (i+1), residue);

        free(B);
        free(Delta);
        free(Y);
        free(R);
    }        
    //pivoteia a matriz inversa para sair corretamente
    pivotMatrix(invMatrix, pRecord);    

    if(config.outputFile == NULL){
        printf("# Tempo LU: %.17g\n", luDecompTime);
        printf("# Tempo iter: %.17g\n", totalLsSolveTime/config.iterationCount);
        printf("# Tempo residuo: %.17g\n", totalResidueTime/config.iterationCount);
        printf("#\n");
        printMatrix(invMatrix, N);        
    }
    else{
        fprintf(config.outputFile, "# Tempo LU: %.17g\n", luDecompTime);
        fprintf(config.outputFile, "# Tempo iter: %.17g\n", totalLsSolveTime/config.iterationCount);
        fprintf(config.outputFile, "# Tempo residuo: %.17g\n", totalResidueTime/config.iterationCount);
        fprintf(config.outputFile, "#\n");
        printMatrixToFile(invMatrix, N, config.outputFile);
    }

    return 0;
}
