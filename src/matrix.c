#include "matrix.h"
#include "util.h"


/**
 * Função pivotMatrix:
 *  Objetivo: Pivotar a matriz de entrada com base em um histórico de pivots
 *  Entrada: - m -> Um vetor que representa uma matriz
 *           - pRecord -> A estrutura de Dados que guarda o historico de pivots realizados
 */
void pivotMatrix(double* m, pivotsRecord* pRecord){
    int i;

    for(i = 0; i < pRecord->count; i++){
        swapRows(m, pRecord->pivots[i]->fromIndex, pRecord->pivots[i]->toIndex, 0, N);
    }
        
}

/**
 * Função swapRows:
 *  Objetivo: Trocar o valor de duas linhas da matriz A a partir da coluna startColumn
 *  Entrada: - A -> Um vetor que representa uma matriz
 *           - fromI -> A linha que vai ser trocada
 *           - toI -> A outra linha que vai ser trocada
 *           - startColumn -> A coluna inicial para a troca
 *           - n -> A ordem da matriz de entrada
 */
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

/**
 * Função pivotRows:
 *  Objetivo: Encontrar o maior elemento em modulo da coluna j e trocar as linhas fazendo com que o maior numero em modulo vire pivo
 *  Entrada: - matrix -> Um vetor que representa uma matriz
 *           - j -> A posicao do pivo a ser encontrado
 *           - n -> A ordem da matriz de entrada
 *           - pRecord -> A estrutura de dados na qual o historico de pivots sera salvado 
 */
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

/**
 * Função copyMatrix:
 *  Objetivo: Copiar os valores da Matriz A para a Matriz B
 *  Entrada: - A -> A matriz a ser copiada
 *           - B -> A matriz a qual vai receber os valores de A
 *           - n -> A ordem das matrizes de entrada 
 */
void copyMatrix(double *A, double *B, unsigned int n)
{
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            getVal(B, i, j) = getVal(A, i, j);
}


/**
 * Função luDecomposition:
 *  Objetivo: Realizar a fatoracao LU da matriz A recebida na entrada
 *  Entrada: - A -> A matriz a ser fatorada
 *           - pRecord -> A estrutura de dados na qual o historico de pivots sera salvado 
 *           - n -> A ordem das matrizes de entrada 
 */
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


/**
 * Função forwardSubstitution:
 *  Objetivo: Solucionar o sistema Ly=B de acordo com a abordagem de substituicao progressiva
 *  Entrada: - L -> A matriz de coeficientes gerada pela fatoracao LU
 *           - B -> Uma matriz contendo em cada coluna o vetor b
 *           - yOrder -> A ordem da matriz B
 *           - n -> A ordem da matriz A 
 */
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

/**
 * Função backwardSubstitution:
 *  Objetivo: Solucionar o sistema Ux=Y de acordo com a abordagem de substituicao regressiva
 *  Entrada: - U -> A matriz triangular superior gerada pela fatoracao LU
 *           - Y -> Uma matriz contendo em cada coluna o vetor y
 *           - xOrder -> A ordem da matriz A
 *           - n -> A ordem da matriz U 
 */
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

/**
 * Função multiplyMatrix:
 *  Objetivo: Realizar a multiplicacao das matrizes A e B e retornar uma nova matriz contendo o resultado
 *  Entrada: - A -> A matriz a ser multiplicada na forma AxB
 *           - B -> A outra matriz a ser multiplicada 
 *           - n -> A ordem das matrizes A e B 
 */
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

/**
 * Função addMatrix:
 *  Objetivo: Realizar a operacao A=A+B
 *  Entrada: - A -> A matriz a ser somada e a que vai receber o resultado
 *           - B -> A outra matriz a ser somada
 *           - n -> A ordem das matrizes A e B 
 */
void AddMatrix(double *A, double *B, unsigned int n){
    int i, j;        

    for(i = 0; i < n; i ++){
        for(j = 0; j < n; j++){
            getVal(A, i, j) += getVal(B, i, j);
        }
    }        
}

/**
 * Função getResidue:
 *  Objetivo: Calcular o erro em cada posicao da matriz e o residuo
 *  Entrada: - Xcandidate -> O candidato a solucao
 *           - X -> A solucao esperada
 *           - residue -> variavel onde sera retornado o residuo
 *           - n -> ordem da matriz de entrada
 */
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