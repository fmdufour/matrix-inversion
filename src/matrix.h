typedef struct LU{
    double* L;
    double* U;
} LU;


void copyMatrix(double *A, double *B, unsigned int n);
void pivot(double* matrix, int j, unsigned int n);
LU* luDecomposition(double* A,unsigned int n);