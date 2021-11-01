#include <fstream>
#include <iostream>
#include "math.h"
#include "common.h"
#define MEM 100000 // memory capacity

int nX, nY; //
int N; // matrix dimension
int i, j, k, l, find, factRet, iter, stop, istep = 0; // counters and flags
int *port, *pos; // portrait and positions of matrix elements
real eps; // accuracy 
real result = 0.0; // dot product result
real *gridX, *gridY; // x-axis and y-axis grids
real *global; // memory location pointer
real *Lo, *Uo; // lower and upper triangles of the original matrix
real *L, *U; // lower and upper triangles of the conditional matrix
real *diOrig, *diCond; // diagonals of original and conditional matrices
real *x, *f; // solution vector and right side vector
real *sout, *s, *r, *p, *z, *q; // auxiliary vectors

int main()
{
    readParams(); // reading parameters
    readGrid(); // reading gridX and gridY
    config(); // configuration of pointers and memory allocation
    //globalAsm(); // assembling the global matrix
    // ...
    factRet = factorization(); // perform LU decomposition
    std::cout << factRet; // удалить
    // здесь проверка корректности факторизации
    //lowerGauss(); // solving the lower triangle
    //upperGauss(); // solving the upper triangle
    //makeResult(); // ...
    return 0;
}



void readParams()
{
    std::ifstream input("parameters.txt");
    input >> nX >> nY;
    N = nX * nY;
    input.close();
}

void readGrid()
{
    std::ifstream inputX("gridX.txt");
    for (i = 0; i < nX; i++)
        inputX >> gridX[i];
    inputX.close();

    std::ifstream inputY("gridY.txt");
    for (i = 0; i < nY; i++)
        inputY >> gridY[i];
    inputY.close();
}

void config()
{
    global = new double[MEM]; // allocation
    memset(global, 0, MEM * sizeof(double)); // nullification

    port = (int*)global;
    pos = (int*)(global + N + 1);

    for (i = 0; i < N + 1; i++) // generating the number of matrix elements
    {
        port[i] = istep;
        istep += i;
    }

    istep = 0;
    for (i = 0; i < N; i++) // generating matrix element positions
        for (k = 0; k < i; k++)
        {
            pos[istep] = k;
            istep++;
        }

    diCond = q + N;
    diOrig = global + 3 * (port[N]) + N + 1;
    L = diCond + N;
    U = L + port[N];
    Lo = global + N + 1 + port[N];
    Uo = global + 2 * (port[N]) + N + 1;
    x = U + port[N];
    f = diOrig + N
    s = x + N;
    sout = s + N;
    r = f + N;
    z = r + N;
    p = z + N;
    q = p + N;
    gridX = sout + N;
    gridY = gridX + nX;
}

int factorization()
{
    for (i = 0; i < N; i++)
    {
        for (j = port[i]; j < port[i + 1]; j++)
        {
            L[j] = (Lo[j] - dotProduct(i, pos[j]));
            U[j] = (Uo[j] - dotProduct(pos[j], i)) / diCond[pos[j]];
        }
        diCond[i] = diOrig[i] - dotProduct(i, i);
    }
    return 0;
}

real dotProduct(int i, int j) // scalar product (for factorization only)
{
    if (i == j)
    {
        for (k = port[i]; k < port[i + 1]; k++)
            result += U[k] * L[k];
    }
    else
    {
        if (i > j) // upper triangle
        {
            for (k = port[j]; k < port[j + 1]; k++)
            {
                find = 0;
                for (l = port[i]; l < port[i + 1] && find == 0; l++)
                {
                    if (pos[l] == pos[k])
                    {
                        result += U[k] * L[l];
                        find = 1;
                    }
                }
            }
        }
        else // lower triangle
        {
            for (l = port[i]; l < port[i + 1]; l++)
            {
                find = 0;
                for (k = port[j]; k < port[j + 1] && find == 0; k++)
                {
                    if (pos[l] == pos[k])
                    {
                        result += U[k] * L[l];
                        find = 1;
                    }
                }
            }
        }

    }
    return result;
}

