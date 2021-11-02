#include <fstream>
#include <iostream>
#include "math.h"
#include "common.h"
#define MEM 100000 // memory capacity

int nX, nY; // number of x-axis and y-axis elements
int N; // matrix dimension
int i, j, k, l, iter, stop, istep = 0; // counters
int factRet; // return value of factorization()
int *port, *pos; // portrait and positions of matrix elements
real eps; // accuracy 
real *gridX, *gridY; // x-axis and y-axis grids
real *global; // memory location pointer
real *Lo, *Uo; // lower and upper triangles of the original matrix
real *L, *U; // lower and upper triangles of the conditional matrix
real *diOrig, *diCond; // diagonals of original and conditional matrices
real *x, *f; // solution vector and right side vector
real *sout, *s, *r, *p, *z, *q; // auxiliary vectors

int main()
{
    config();
    //globalAsm();
    factRet = factorization();
    lowerGauss(f, r);
    upperGauss(r, z);
    multMatVect(z, q);
    lowerGauss(q, p);
    //makeResult();
    return 0;
}

void config() // configuration of pointers and memory allocation
{
    global = new double[MEM]; // allocation
    memset(global, 0, MEM * sizeof(double)); // nullification
    readSize();

    port = (int*)global;
    pos = (int*)(global + N + 1);

    for (i = 0; i < N + 1; i++) // generating the number of matrix elements
    {
        port[i] = istep;
        istep += i;
    }

    istep = 0;
    for (i = 0; i < N; i++) // generating matrix element positions
    {
        for (k = 0; k < i; k++)
        {
            pos[istep] = k;
            istep++;
        }
    }

    Lo = global + N + 1 + port[N];
    Uo = global + 2 * (port[N]) + N + 1;
    diOrig = global + 3 * (port[N]) + N + 1;
    f = diOrig + N;
    r = f + N;
    z = r + N;
    p = z + N;
    q = p + N;
    diCond = q + N;
    L = diCond + N;
    U = L + port[N];
    x = U + port[N];
    s = x + N;
    sout = s + N;
    gridX = sout + N;
    gridY = gridX + nX;

    readGrid();
}

//void globalAsm() // assembling the global matrix
//{
//
//}

int factorization() // perform LU decomposition
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

void addGlobal(int i, int j, real elem) // adds elements to the global matrix
{
    if (i == j)
        diOrig[i] += elem;

    else if (i > j)
    {
        for (k = port[i]; k < port[i + 1]; k++)
            if (pos[k] == j)
                Lo[k] += elem;
    }

    else
    {
        for (k = port[j]; k < port[j + 1]; k++)
            if (pos[k] == i)
                Uo[k] += elem;
    }
}