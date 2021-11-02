#include <fstream>
#include "common.h"
#define MEM 100000

int nX, nY; // number of x-axis and y-axis elements
int N; // matrix dimension
int factRet; // return value of factorization()
int *port, *pos; // portrait and positions of matrix elements
real eps; // accuracy 
real *gridX, *gridY; // x-axis and y-axis grids
real *BUF; // memory location pointer
real *Lo, *Uo; // lower and upper triangles of the original matrix
real *L, *U; // lower and upper triangles of the conditional matrix
real *diOrig, *diCond; // diagonals of original and conditional matrices
real *x, *f; // solution vector and right side vector
real *sout, *s, *r, *p, *z, *q; // auxiliary vectors

int main()
{
    BUF = new double[MEM]; // dynamic memory allocation
    memset(BUF, 0, MEM * sizeof(double)); // nullification of the buffer

    readSize();
    config();
    readGrid();
    //globalAsm();
    //boundaryCond();
    factRet = factorization();
    lowerGauss(f, r);
    upperGauss(r, z);
    multMatVect(z, q);
    lowerGauss(q, p);
    //makeResult();
    return 0;
}

void config() // configuration of pointers
{
    port = (int*)BUF;
    pos = (int*)(BUF + N + 1);

    int i, k, step = 0;
    for (i = 0; i < N + 1; i++) // generating the number of matrix elements
    {
        port[i] = step;
        step += i;
    }

    step = 0;
    for (i = 0; i < N; i++) // generating matrix element positions
    {
        for (k = 0; k < i; k++)
        {
            pos[step] = k;
            step++;
        }
    }

    Lo = BUF + N + 1 + port[N];
    Uo = BUF + 2 * (port[N]) + N + 1;
    diOrig = BUF + 3 * (port[N]) + N + 1;
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
}

//void globalAsm() // assembling the global matrix
//{
//
//}

//void boundaryCond() // working with boundary conditions
//{
//
//}

void addGlobal(int i, int j, real elem) // adding elements to the global matrix
{
    int k;

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