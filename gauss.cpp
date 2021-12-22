#include "common.h"

void forwardGauss(double *in, double *out) // solving the lower triangle
{
    int i, j;
    double result;

    for (i = 0; i < N; i++)
    {
        result = 0;
        for (j = ig[i]; j < ig[i + 1]; j++)
            result += L[j] * out[jg[j]];

        out[i] = (in[i] - result) / di[i];
    }
}

void backwardGauss(double *in, double *out) // solving the upper triangle
{
    int i, j;

    for (i = 0; i < N; i++)
        out[i] = in[i];

    for (i = N - 1; i >= 0; i--)
        for (j = ig[i]; j < ig[i + 1]; j++)
            out[jg[j]] -= U[j] * out[i];
}

void factorization() // perform LU decomposition
{
    int i, j;

    for (i = 0; i < N; i++)
    {
        for (j = ig[i]; j < ig[i + 1]; j++)
        {
            L[j] = (ggl[j] - sum(i, jg[j]));
            U[j] = (ggu[j] - sum(jg[j], i)) / di[jg[j]];
        }
        di[i] = diOrig[i] - sum(i, i);
    }
}