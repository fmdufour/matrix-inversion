#include "matrix.h"
#include "util.h"

void pivot(double* matrix, int j, unsigned int n){
    int i = j;
    double temp;    
    double max = abs(getVal(matrix, i, j));
    int maxIndex = i;
    
    for(i++; i < n; i++){
        if(abs(getVal(matrix, i, j)) > max){
            max = abs(getVal(matrix, i, j));
            maxIndex = i;
        }
    }

    if(maxIndex != j){        
        for(i = j; i < N; i++){
            temp = getVal(matrix, j, i);
            getVal(matrix, j, i) = getVal(matrix, maxIndex, i);
            getVal(matrix, maxIndex, i) = temp;
        }
    }    
}

void copyMatrix(double *A, double *B, unsigned int n){
    int i,j;
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            getVal(B, i, j) = getVal(A, i, j);
}

LU* luDecomposition(double* A,unsigned int n){
    int i, j, colIndex;
    LU* lu = malloc(sizeof(LU));

    lu->L = allocateMatrix(n);
    lu->U = allocateMatrix(n);

    copyMatrix(A, lu->U, n);    
    
    for(j = 0; j < N; j++){
        pivot(lu->U, j, N);
        for(i = j+1; i < N; i++){
            getVal(lu->L, i, j) = getVal(lu->U, i, j)/getVal(lu->U, j, j);
            
            for(colIndex = j; colIndex < n; colIndex++){
                getVal(lu->U, i, colIndex) 
                    = getVal(lu->U, i, colIndex) - getVal(lu->L, i, j) * getVal(lu->U, j, colIndex);
            }
        }
    }
    
    return lu;
}

void forwardSubstitution(double* L, double* Y, int colIndex, unsigned n){
    int i, k;
    double sum;
    double* B = getIdentityColumn(colIndex, n);
    
    getVal(Y, 0, colIndex) = B[0]/getVal(L, 0, 0);      

    for(i = 1; i < n; i++){
        sum = 0;
        for(k = 0; k < (i-1); k++)
            sum += sum + getVal(L, i, k) * getVal(Y, k, colIndex);
        
        getVal(Y, i, colIndex) = (B[i] - sum) / getVal(L, i, i);
    }    
}

void backwardSubstitution(double* U, double* X, double* Y, int colIndex, unsigned n){
    int i, k;
    double sum;    
    
    getVal(X, n-1, colIndex) = getVal(Y, n-1, colIndex)/ getVal(U, n-1, n-1);
    //getVal(Y, n-1, colIndex) = B[0]/getVal(L, 0, 0);      

    for(i = n-2; i >= 0; i--){
        sum = 0;
        for(k = i + 1; k < n; k++)
            sum += sum + getVal(U, i, k) * getVal(X, k, colIndex);
        
        getVal(X, i, colIndex) = (getVal(Y, i, colIndex) - sum) / getVal(U, i, i);
    }
    // xn,j = yn,j / un,n;
    // for i := n-1 downto 1 do /* process elements of that column */
    
    //     begin
    //         sum := 0; /* solving for xi on the current column */
    //         for k := i+1 to n do
    //             sum := sum + ui,k× xk,j;
    //         xi,j = (yi,j - sum)/ui,i; end 
}




















