#include "util.h"

int main(int argc, char* argv[]){    
    int i;
    configuration config;    
    int receivedIterationFlag = 0;

    for(i = 0; i < argc; i++){
        if(strcmp(argv[i], "-e") == 0){            
            if(++i == argc){
                printErrorExit("Informe o nome do Arquivo de entrada\n-o <nomeArquivoEntrada>\n");
            }            
            config.inputFile = argv[i];
        }
        else if(strcmp(argv[i], "-o") == 0){            
            if(++i == argc){
                printErrorExit("Informe o nome do Arquivo de saída\n-o <nomeArquivoSaida>\n");
            }            
            config.outputFile = argv[i];
        }        
        else if(strcmp(argv[i], "-i") == 0){            
            if(++i == argc){
                printErrorExit("Informe o número de Iteracoes do refinamento\n-i <numIteracoes>\n");
            }
            receivedIterationFlag = 1;
            config.iterationCount = atoi(argv[i]);
        }
        else if(strcmp(argv[i], "-r") == 0){
            if(++i == argc){
                printErrorExit("Informe o tamanho da matriz aleatoria\n-r <tamMatrizAleatoria>\n");
            }            
            config.randomMatrixSize = atoi(argv[i]);
        }            
    }

    if(!receivedIterationFlag){
        printErrorExit("Informe o número de Iteracoes do refinamento\n-i <numIteracoes>\n");
    }

    printf("Configurações recebidas:\n");
    printf("Arquivo de Entrada: %s\n", config.inputFile);
    printf("Arquivo de Saida: %s\n", config.outputFile);
    printf("Tam Matriz aleatoria: %d\n", config.randomMatrixSize);
    printf("Numero de Iteracoes Refinamento: %d\n", config.iterationCount);

    srand(20172);
    double** matrix = generateSquareRandomMatrix(10);
    printMatrix(matrix, 10);    
    return 0;
}