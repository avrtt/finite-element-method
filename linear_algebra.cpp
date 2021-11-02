#include "common.h"

real dotProduct(int i, int j) // scalar product (for factorization)
{
    int k, l, find;
    real result = 0.0;

    if (i == j)
    {
        for (k = port[i]; k < port[i + 1]; k++)
            result += U[k] * L[k];
    }

    else if (i > j) // upper triangle
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
    return result;
}

real multVect(real *vect1, real *vect2) // scalar multiplication of two vectors
{
    int i;
    real result = 0;

    for (i = 0; i < N; i++)
    {
        result += vect1[i] * vect2[i];
    }
    return result;
}

void multMatVect(real *in, real *out) // matrix-vector product
{
    int i, j;
    real *outNew = new real[N];

    for (i = 0; i < N; i++)
    {
        outNew[i] = diOrig[i] * in[i];
        for (j = port[i]; j < port[i + 1]; j++)
        {
            outNew[i] += Lo[j] * in[pos[j]];
            outNew[pos[j]] += Uo[j] * in[i];
        }
    }

    for (i = 0; i < N; i++)
        out[i] = outNew[i];

    delete[] outNew;
}