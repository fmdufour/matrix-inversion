#include <math.h>


/**
 * As matrizes sao armazenas em vetores, para otimizar o acesso fazendo com que seja armazenado em areas proximas da memoria
 * Estrutura de dados para armazenar as matrizes L e U
 */
typedef struct LU{
    double* L;
    double* U;
} LU;


/**
 * Estrutura de dados para armazenar o pivot, guarda o i e o j das linhas trocadas
 */
typedef struct pivot{
    int fromIndex;
    int toIndex;
} pivot;


/**
 * Estrutura de dados para armazenar o historico de pivots realizados, guarda o i e o j das linhas trocadas
 */
typedef struct pivotsRecord{
    pivot** pivots;
    int count;
} pivotsRecord;

void AddMatrix(double *A, double *B, unsigned int n);
void pivotMatrix(double* m, pivotsRecord* pRecord);
double* getResidue(double* B, double* idMatrix, double *residue, unsigned int n);
void swapRows(double* A, int fromI, int toI, int startColumn, int n);
double* multiplyMatrix(double* A, double* B, unsigned int n);
void copyMatrix(double *A, double *B, unsigned int n);
void pivotRows(double* matrix, int j, unsigned int n, pivotsRecord* pivots);
LU* luDecomposition(double* A, pivotsRecord* pRecord, unsigned int n);
double getDeterminant(double* m, unsigned int n);