#include <fstream>
#include "math.h"
#include "common.h"
#define MEM 100000 // memory capacity

int nX, nY;
int N; // matrix dimension
int i, k, istep = 0; // counters
int *ig, *jg; // portrait and positions of matrix elements
real eps; // accuracy 
real *gridX, *gridY; // x-axis and y-axis grids
real *global; // memory location pointer
real *ggl, *ggu; // lower and upper triangles of the original matrix
real *L, *U; // lower and upper triangles of the conditional matrix
real *di, *diag; // diagonals of original and conditional matrices
real *x, *f; // solution vector and right side vector
real *r, *p, *z, *q, *s, *sout; // auxiliary vectors

int main()
{
    readParams(); // reading parameters
    readGrid(); // reading gridX and gridY
    config(); // configuration of pointers and memory allocation
    //globalAsm(); // assembling the global matrix
    // ...
    //factorization(); // perform LU decomposition
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

    ig = (int*)global;
    jg = (int*)(global + N + 1);

    for (i = 0; i < N + 1; i++) // generating the number of matrix elements
    {
        ig[i] = istep;
        istep += i;
    }

    istep = 0;
    for (i = 0; i < N; i++) // generating matrix element positions
        for (k = 0; k < i; k++)
        {
            jg[istep] = k;
            istep++;
        }

    ggl = global + N + 1 + ig[N];
    ggu = global + 2 * (ig[N]) + N + 1;
    di = global + 3 * (ig[N]) + N + 1;
    f = di + N;
    r = f + N;
    z = r + N;
    p = z + N;
    q = p + N;
    diag = q + N;
    L = diag + N;
    U = L + ig[N];
    x = U + ig[N];
    s = x + N;
    sout = s + N;
    gridX = sout + N;
    gridY = gridX + nX;
}

