#include "common.h"

void multVectByA(double *in, double *out) // multiplication of a given vector by the matrix A 
{
    int i, j;

    for (i = 0; i < N; i++)
    {
        out[i] = diOrig[i] * in[i];
        for (j = ig[i]; j < ig[i + 1]; j++)
        {
            out[i] += ggl[j] * in[jg[j]];
            out[jg[j]] += ggu[j] * in[i];
        }
    }
}

double multVectByVect(double *vect1, double *vect2) // scalar multiplication of two vectors
{
    int i;
    double result = 0;

    for (i = 0; i < N; i++)
        result += vect1[i] * vect2[i];

    return result;
}

double sum(int i, int j) // sum of L and U products (for factorization)
{
    int k, l, flag;
    double result = 0.0;

    if (i > j) 
    {
        for (k = ig[j]; k < ig[j + 1]; k++)
        {
            flag = 0;
            for (l = ig[i]; l < ig[i + 1] && flag == 0; l++)
            {
                if (jg[l] == jg[k])
                {
                    result += U[k] * L[l];
                    flag = 1;
                }
            }
        }
    }

    else if (i < j)
    {
        for (l = ig[i]; l < ig[i + 1]; l++)
        {
            flag = 0;
            for (k = ig[j]; k < ig[j + 1] && flag == 0; k++)
            {
                if (jg[l] == jg[k])
                {
                    result += U[k] * L[l];
                    flag = 1;
                }
            }
        }
    }

    else 
        for (k = ig[i]; k < ig[i + 1]; k++)
            result += U[k] * L[k];

    return result;
}

