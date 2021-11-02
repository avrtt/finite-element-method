#include "math.h"
#include "common.h"

// this section allows you to specify parameter functions depending on x and y

real analytical(real x, real y)
{
    return 0;
}

real getF(real x, real y)
{
    return x * x - y * y;
}

real getLambda(real x, real y)
{
    return 1.0;
}

real getGamma(real x, real y)
{
    return 0;
}