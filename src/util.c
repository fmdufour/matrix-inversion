#include "util.h"


/**
 * Função getIdentityMatrix:
 *  Objetivo: Retornar uma matriz identidade de ordem n
 *  Entrada: - n -> A ordem da matriz a ser gerada
 *  Saida: A matriz identidade 
 */
double* getIdentityMatrix(unsigned int n){
    int i,j;
    
    double* identity = allocateMatrix(n);    

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            getVal(identity, i, j) = (i == j ? 1.0 : 0.0);                        
        }        
    }
        
    return identity;
}

/**
 * Função allocateMatrix:
 *  Objetivo: Alocar um vetor para representar uma matriz de ordem n
 *  Entrada: - n -> ordem da matriz a ser alocada
 *  Saida: O vetor alocado para a matriz
 */
double* allocateMatrix(unsigned int n){
    double *mat = NULL;
    /* return NULL if memory allocation fails */
    if (!(mat = (double *)malloc(n * n * sizeof(double))))
        return (NULL);

    return mat;
}

/**
 * Função generateSquareRandomMatrix:
 *  Objetivo: Gerar um vetor para representar uma matriz de ordem n com valores aleatorios
 *  Entrada: - n -> ordem da matriz a ser alocada
 *  Saida: O vetor alocado para a matriz com os valores aleatorios
 */
double* generateSquareRandomMatrix(unsigned int n)
{
    double *mat = NULL;

    mat = allocateMatrix(n);

    /* generate a randomly initialized matrix in row-major order */
    double *ptr = mat;
    double *end = mat + n * n;
    double invRandMax = 1.0 / (double)RAND_MAX;

    while (ptr != end){
        *ptr++ = (double)rand() * invRandMax;
    }

    return mat;
}

/**
 * Função readMatrixFromStdIn:
 *  Objetivo: retornar um vetor para representar uma matriz lida do stdin 
 *  Saida: O vetor alocado para a matriz
 */
double* readMatrixFromStdIn(){
    double *mat = NULL;
    int i;

    scanf("%d", &N);    

    mat = allocateMatrix(N);    
    
    for(i = 0; i < N*N; i++)
        scanf("%lf", &mat[i]);

    return mat;
}

/**
 * Função readMatrixFromFile:
 *  Objetivo: retornar um vetor para representar uma matriz lida de um arquivo
 *  Entrada: Estrutura de dados de configuracao a qual contem o FILE*
 *  Saida: O vetor alocado para a matriz 
 */
double* readMatrixFromFile(configuration config){
    double *mat = NULL;
    int i;

    fscanf(config.inputFile, "%d", &N);    

    mat = allocateMatrix(N);    
    
    for(i = 0; i < N*N; i++)
        fscanf(config.inputFile, "%lf", &mat[i]);

    return mat;
}

/**
 * Função timestamp:
 *  Objetivo: retorna um timestamp do momento 
 *  Saida: O timestamp 
 */
double timestamp(void)
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double)(tp.tv_sec * 1000.0 + tp.tv_usec / 1000.0));
}


/**
 * Função printMatrix:
 *  Objetivo: Imprimir no stdout uma matriz 
 */
void printMatrix(double *matrix, unsigned int n)
{
    int i, j;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            printf("%.17g ", getVal(matrix, i, j));
        }
        printf("\n");
    }
}

/**
 * Função printMatrix:
 *  Objetivo: Imprimir uma matriz de ordem n para o arquivo 
 */
void printMatrixToFile(double* matrix, unsigned int n, FILE* file)
{
    int i, j;

    fprintf(file, "%d\n", N);
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            fprintf(file, "%.17g ", getVal(matrix, i, j));
        }
        fprintf(file, "\n");
    }
}

/**
 * Função printErrorExit:
 *  Objetivo: imprimir no stderr uma mensagem e encerra o programa com codigo 1
 */
void printErrorExit(char *msg)
{
    fprintf(stderr, "%s", msg);
    exit(1);
}

/**
 * Função readConfiguration:
 *  Objetivo: alocar a estrutura de dados configuration com base nos parametros recebidos por linha de comando
 * verificando tambem as regras para os mesmos
 */
configuration readConfiguration(int argc, char *argv[])
{
    int i;
    configuration config;
    config.randomMatrixSize = 0;
    config.inputFile = NULL;
    config.outputFile = NULL;
    config.inputFilePath = NULL;
    config.outputFilePath = NULL;
    int receivedIterationFlag = 0;

    for (i = 0; i < argc; i++)
    {
        if (strcmp(argv[i], "-e") == 0)
        {
            if (++i == argc)
            {
                printErrorExit("Informe o nome do Arquivo de entrada\n-o <nomeArquivoEntrada>\n");
            }
            config.inputFilePath = argv[i];
            config.inputFile = fopen(config.inputFilePath, "r");

            if (config.inputFile == NULL)
                printErrorExit("Arquivo de entrada invalido\n");
        }
        else if (strcmp(argv[i], "-o") == 0)
        {
            if (++i == argc)
            {
                printErrorExit("Informe o nome do Arquivo de saída\n-o <nomeArquivoSaida>\n");
            }
            config.outputFilePath = argv[i];
            config.outputFile = fopen(config.outputFilePath, "w");

            if (config.outputFile == NULL)
                printErrorExit("Arquivo de saida invalido\n");
        }
        else if (strcmp(argv[i], "-i") == 0)
        {
            if (++i == argc)
            {
                printErrorExit("Informe o número de Iteracoes do refinamento\n-i <numIteracoes>\n");
            }
            receivedIterationFlag = 1;
            config.iterationCount = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "-r") == 0)
        {
            if (++i == argc)
            {
                printErrorExit("Informe o tamanho da matriz aleatoria\n-r <tamMatrizAleatoria>\n");
            }
            config.randomMatrixSize = atoi(argv[i]);

            if(config.randomMatrixSize == 0){
                printErrorExit("Informe o tamanho da matriz aleatoria\n-r <tamMatrizAleatoria>\n");
            }

            N = config.randomMatrixSize;
        }
    }

    if (!receivedIterationFlag){
        printErrorExit("Informe o número de Iteracoes do refinamento\n-i <numIteracoes>\n");
    }

    if(config.inputFile == NULL && config.randomMatrixSize == 0){
        printErrorExit("Informe o tamanho da matriz aleatoria\n-r <tamMatrizAleatoria>\n");
    }

    return config;
}

/**
 * Função diffInSeconds:
 *  Objetivo: Calcula a diferenca entre o momento e o timestamp recebido por pametro em segundos
 */
double diffInSeconds(double start){
    return (timestamp() - start)/1000;
}
