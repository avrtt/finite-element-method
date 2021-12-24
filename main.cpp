#include "common.h"

int nX, nY; // number of x-axis and y-axis elements
int N; // linear system dimension
int i, k, n, m, step = 0, flag = 0; // counters and flags
int *ig, *jg; // line beginning indexes and positions of matrix elements
double *gridX, *gridY; // x-axis and y-axis grids
double *edges; // types of boundary conditions on four sides
double *ggl, *ggu; // lower and upper triangles of the global matrix
double *L, *U; // lower and upper triangles of the factorized global matrix
double *diOrig, *di; // diagonals of the global matrix before and after factorization
double *x, *f; // solution vector and right side vector
double *r, *p, *z, *q; // auxiliary vectors for LOS
double alpha, beta, breakValue; // LOS coefficients and exit value
double xPoint, yPoint, xUpperPoint, yUpperPoint, hx, hy; // finite element parameters
double lambda, gamma; // area parameters
double temp, hx_2, hy_2, el1, el2, el3, el4; // variables for elements of G and C

int main()
{
    readSize(); 
    config();
    readGrid();
    globalAsm();
    boundaryCond();
    factorization();
    los();
    printResult();

    return 0;
}

void globalAsm() // building local matrices and assembling the global matrix
{
    int node[4]; // node numbers
    double G[4][4], C[4][4]; // stiffness and mass (with gamma = 1) matrices
    double fl[4], b[4]; // right side values and local vector
    double A[4][4]; // local matrix 

    for (k = 0; k < nY - 1; k++)
        for (i = 0; i < nX - 1; i++)
        {
            // points defining the finite element
            xPoint = gridX[i];
            xUpperPoint = gridX[i + 1];
            yPoint = gridY[k];
            yUpperPoint = gridY[k + 1];

            // horizontal and vertical steps
            hx = xUpperPoint - xPoint;
            hy = yUpperPoint - yPoint;

            lambda = getLambda(xPoint + hx / 2.0, yPoint + hy / 2.0);
            gamma = getGamma(xPoint + hx / 2.0, yPoint + hy / 2.0);

            hx_2 = hx * hx;
            hy_2 = hy * hy;
            temp = 1.0 / (hx * hy);
            el1 = (hx_2 + hy_2) * temp / 3;
            el2 = (hx_2 - 2 * hy_2) * temp / 6;
            el3 = - (2 * hx_2 - hy_2) * temp / 6;
            el4 = - (hx_2 + hy_2) * temp / 6;

            G[0][0] = el1, G[0][1] = el2, G[0][2] = el3, G[0][3] = el4;
            G[1][0] = el2, G[1][1] = el1, G[1][2] = el4, G[1][3] = el3;
            G[2][0] = el3, G[2][1] = el4, G[2][2] = el1, G[2][3] = el2;
            G[3][0] = el4, G[3][1] = el3, G[3][2] = el2, G[3][3] = el1; 

            /*
            
            lambda1 = getLambda(xPoint, yPoint);
            lambda2 = getLambda(xPoint + hx / 2.0, yPoint);
            lambda3 = getLambda(xPoint + hx, yPoint);
            lambda4 = getLambda(xPoint, yPoint + hy / 2.0);
            lambda5 = getLambda(xPoint + hx / 2.0, yPoint + hy / 2.0);
            lambda6 = getLambda(xPoint + hx, yPoint + hy / 2.0);
            lambda7 = getLambda(xPoint, yPoint + hy);
            lambda8 = getLambda(xPoint + hx / 2.0, yPoint + hy);
            lambda9 = getLambda(xPoint + hx, yPoint + hy);
            gamma = getGamma(xPoint + hx / 2.0, yPoint + hy / 2.0);

            G[0][0] = lambda1 * ... + ... + lambda9 * ...;
            G[1][0] = lambda1 * ... + ... +  lambda9 * ...;
            ... 
            
            */

            el1 = hx * hy / 9.0;
            el2 = hx * hy / 18.0;
            el3 = el2;
            el4 = hx * hy / 36.0;

            C[0][0] = el1, C[0][1] = el2, C[0][2] = el3, C[0][3] = el4;
            C[1][0] = el2, C[1][1] = el1, C[1][2] = el4, C[1][3] = el3;
            C[2][0] = el3, C[2][1] = el4, C[2][2] = el1, C[2][3] = el2;
            C[3][0] = el4, C[3][1] = el3, C[3][2] = el2, C[3][3] = el1;

            fl[0] = getF(xPoint, yPoint);
            fl[1] = getF(xUpperPoint, yPoint);
            fl[2] = getF(xPoint, yUpperPoint);
            fl[3] = getF(xUpperPoint, yUpperPoint);

            b[0] = C[0][0] * fl[0] + C[0][1] * fl[1] + C[0][2] * fl[2] + C[0][3] * fl[3];
            b[1] = C[1][0] * fl[0] + C[1][1] * fl[1] + C[1][2] * fl[2] + C[1][3] * fl[3];
            b[2] = C[2][0] * fl[0] + C[2][1] * fl[1] + C[2][2] * fl[2] + C[2][3] * fl[3];
            b[3] = C[3][0] * fl[0] + C[3][1] * fl[1] + C[3][2] * fl[2] + C[3][3] * fl[3];

            node[0] = nX * k + i;
            node[1] = nX * k + i + 1;
            node[2] = nX * (k + 1) + i;
            node[3] = nX * (k + 1) + i + 1;

            // adding elements to the global matrix
            for (n = 0; n < 4; n++)
            {
                f[node[n]] += b[n]; // assembling the global vector

                for (m = 0; m < 4; m++)
                {
                    A[n][m] = lambda * G[n][m] + gamma * C[n][m]; // building a local matrix
                    addGlobal(A[n][m], node[n], node[m]); // adding to the global matrix
                }
            }
        }
}

void addGlobal(double elem, int i, int j)
{
    int k;

    if (i > j)
    {
        for (k = ig[i]; k < ig[i + 1]; k++)
            if (jg[k] == j)
                ggl[k] += elem;
    }

    else if (i < j)
    {
        for (k = ig[j]; k < ig[j + 1]; k++)
            if (jg[k] == i)
                ggu[k] += elem;
    }

    else diOrig[i] += elem;
}

void los()
{
    for (i = 0; i < N; i++)
    {
        x[i] = 0; // solution vector
        r[i] = f[i]; // r = f - Ax
    }

    forwardGauss(r, r); // r = L'r
    backwardGauss(r, z); // z = U'r
    multVectByA(z, p); // p = Az
    forwardGauss(p, p); // p = L'Az

    for (k = 0; k < numIter && flag == 0; k++)
    {
        alpha = multVectByVect(p, r) / multVectByVect(p, p);

        for (i = 0; i < N; i++)
        {
            x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
        }

        backwardGauss(r, q); // q = U'r
        multVectByA(q, q); // q = AU'r
        forwardGauss(q, q); // q = L'AU'r

        beta = - (multVectByVect(p, q) / multVectByVect(p, p));

        for (i = 0; i < N; i++)
        {
            z[i] = beta * z[i] + q[i];
            p[i] = beta * p[i] + q[i];
        }

        breakValue = multVectByVect(r, r);

        if (breakValue < eps)
            flag = 1;
    }
}