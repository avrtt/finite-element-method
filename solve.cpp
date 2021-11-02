#include "common.h"

int factorization() // perform LU decomposition
{
    int i, j;

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

void lowerGauss(real *in, real *out) // solving the lower triangle
{
    int i, j;
    real result;

    for (i = 0; i < N; i++)
    {
        result = 0;
        for (j = port[i]; j < port[i + 1]; j++)
        {
            result += L[j] * out[pos[j]];
        }
        out[i] = (in[i] - result) / diCond[i];
    }
}

void upperGauss(real *in, real *out) // solving the upper triangle
{
    int i, j;

    for (i = 0; i < N; i++)
        out[i] = in[i];

        for (i = N - 1; i >= 0; i--)
        {
            for (j = port[i]; j < port[i + 1]; j++)
                out[pos[j]] -= U[j] * out[i];
        }
}