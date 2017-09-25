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


void swapRows(double* A, int fromI, int toI, int startColumn, int n);
double* multiplyMatrix(double* A, double* B, unsigned int n);
void copyMatrix(double *A, double *B, unsigned int n);
void pivotRows(double* matrix, int j, unsigned int n, pivotsRecord* pivots);
LU* luDecomposition(double* A, pivotsRecord* pRecord, unsigned int n);