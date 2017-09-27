#include <math.h>

typedef struct LU{
    double* L;
    double* U;
} LU;

typedef struct pivot{
    int fromIndex;
    int toIndex;
} pivot;

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