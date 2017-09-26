#include "matrix.h"
#include "util.h"


void pivotMatrix(double* m, pivotsRecord* pRecord){
    int i;

    for(i = 0; i < pRecord->count; i++){
        swapRows(m, pRecord->pivots[i]->fromIndex, pRecord->pivots[i]->toIndex, 0, N);
    }
        
}

void swapRows(double* A, int fromI, int toI, int startColumn, int n){
    int i;
    double temp;
    for (i = startColumn; i < n; i++)
    {        
        temp = getVal(A, fromI, i);
        getVal(A, fromI, i) = getVal(A, toI, i);
        getVal(A, toI, i) = temp;
    }
}

void pivotRows(double *matrix, int j, unsigned int n, pivotsRecord* pRecord)
{
    int i = j;
    double max = fabs(getVal(matrix, i, j));
    int maxIndex = i;
    pivot* p;    

    for (i++; i < n; i++)
    {
        if (fabs(getVal(matrix, i, j)) > max)
        {
            max = fabs(getVal(matrix, i, j));
            maxIndex = i;
        }
    }

    if (maxIndex != j)
    {
        p = malloc(sizeof(pivot));
        p->fromIndex = j;
        p->toIndex = maxIndex;
        pRecord->pivots[pRecord->count] = p;       
        pRecord->count++;

        swapRows(matrix, j, maxIndex, j, n);        
    }
}

void copyMatrix(double *A, double *B, unsigned int n)
{
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            getVal(B, i, j) = getVal(A, i, j);
}

LU *luDecomposition(double *A, pivotsRecord* pRecord, unsigned int n)
{
    int i, j, colIndex;
    LU *lu = malloc(sizeof(LU));

    lu->L = allocateMatrix(n);
    lu->U = allocateMatrix(n);

    for(i = 0; i < N; i ++)
        getVal(lu->L, i, i) = 1.0;

    copyMatrix(A, lu->U, n);

    for (j = 0; j < N; j++)
    {
        pivotRows(lu->U, j, N, pRecord);  
        for (i = j + 1; i < N; i++)
        {
            getVal(lu->L, i, j) = getVal(lu->U, i, j) / getVal(lu->U, j, j);

            for (colIndex = j; colIndex < n; colIndex++)
            {
                getVal(lu->U, i, colIndex) = getVal(lu->U, i, colIndex) - getVal(lu->L, i, j) * getVal(lu->U, j, colIndex);
            }
        }
    }

    return lu;
}

double* forwardSubstitution(double *L, double *B, unsigned int yOrder, unsigned int n)
{
    int i, k, colIndex;
    double sum;    
    double *Y = allocateMatrix(yOrder);

    for (colIndex = 0; colIndex < yOrder; colIndex++)
    {        
        getVal(Y, 0, colIndex) = getVal(B, 0, colIndex) / getVal(L, 0, 0);

        for (i = 1; i < n; i++)
        {
            sum = getVal(B, i, colIndex);
            for (k = 0; k < i; k++)
            {
                sum -= getVal(L, i, k) * getVal(Y, k, colIndex);
            }

            getVal(Y, i, colIndex) = (sum);
        }        
    }
    
    return Y;
}

double* backwardSubstitution(double *U, double *Y, unsigned int xOrder, unsigned int n)
{
    int i, k, colIndex;
    double sum;
    double *X = allocateMatrix(xOrder);

    for (colIndex = 0; colIndex < xOrder; colIndex++)
    {
        getVal(X, n - 1, colIndex) = getVal(Y, n - 1, colIndex) / getVal(U, n - 1, n - 1);
    
        for (i = n - 2; i >= 0; i--)
        {
            sum = getVal(Y, i, colIndex);
            for (k = i + 1; k < n; k++)
                sum -= getVal(U, i, k) * getVal(X, k, colIndex);
    
            getVal(X, i, colIndex) = (sum) / getVal(U, i, i);
        }
    }

    return X;
}

double* multiplyMatrix(double* A, double* B, unsigned int n){
    int i, j, k;
    double sum;
    double* R = allocateMatrix(n);

    for(i = 0; i < n; i ++){
        for(j = 0; j < n; j++){
            sum = 0;
            for(k = 0; k < n; k++){
                sum += getVal(A, i, k) * getVal(B, k, j);
            }
            getVal(R, i, j) = sum;
        }
    }
        
    return R;
}

void AddMatrix(double *A, double *B, unsigned int n){
    int i, j;        

    for(i = 0; i < n; i ++){
        for(j = 0; j < n; j++){
            getVal(A, i, j) += getVal(B, i, j);
        }
    }        
}

double* getResidue(double *Xcandidate, double *X, double* residue, unsigned int n){
    int i,j;    
    double* R;
    *residue = 0;
    R = allocateMatrix(n);
    
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            getVal(R, i , j) = getVal(X, i, j) - getVal(Xcandidate, i, j);
            *residue += pow(getVal(R, i, j), 2);
        }
    }          
    *residue = sqrt(*residue);

    return R;
}