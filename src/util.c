#include "util.h"

double* getIdentityColumn(unsigned int iCol, unsigned int n){
    double* col = malloc(sizeof(double) * n);
    int i;

    for(i = 0; i < n; i++){
        col[i] = iCol == i ? 1.0 : 0.0;        
    }
        
    return col;
}

double* allocateMatrix(unsigned int n){
    double *mat = NULL;
    /* return NULL if memory allocation fails */
    if (!(mat = (double *)malloc(n * n * sizeof(double))))
        return (NULL);

    return mat;
}

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

double* readMatrixFromFile(configuration config){
    double *mat = NULL;
    int i;

    fscanf(config.inputFile, "%d", &N);    

    mat = allocateMatrix(N);    
    
    for(i = 0; i < N*N; i++)
        fscanf(config.inputFile, "%lf", &mat[i]);

    return mat;
}

double timestamp(void)
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double)(tp.tv_sec * 1000.0 + tp.tv_usec / 1000.0));
}

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

void printErrorExit(char *msg)
{
    fprintf(stderr, msg);
    exit(1);
}

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
