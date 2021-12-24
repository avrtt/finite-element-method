#include "math.h"
#include "common.h"

#define MEM 1000000 // memory capacity
int numIter = 1000; // maximum number of LOS iterations 
double eps = 1e-30; // solution accuracy of LOS

double getAnalytical(double x, double y)
{
    return pow(x, 4) + pow(y, 4)
;}

double getF(double x, double y)
{
    return - 12 * x * x - 12 * y * y 
;}

double getLambda(double x, double y)
{
    return 1
;}

double getGamma(double x, double y)
{
    return 0
;}

double getBeta()
{
    return 1
;}


void buildPortrait() // generating the matrix portrait
{
    for (i = 0; i < N + 1; i++)
    {
        ig[i] = step;
        step += i;
    }

    step = 0;
    for (i = 0; i < N; i++)
        for (k = 0; k < i; k++)
        {
            jg[step] = k;
            step++;
        }
}

void config()
{
    double *BUF = new double[MEM]; // dynamic memory allocation
    memset(BUF, 0, MEM * sizeof(double)); // nullification of the buffer

    ig = (int*)(BUF);
    jg = (int*)(BUF + N + 1);

    buildPortrait();

    ggl = BUF + N + ig[N];
    ggu = ggl + N + ig[N];
    diOrig = ggu + N + ig[N];

    r = diOrig + N;
    z = r + N;
    p = z + N;
    q = p + N;

    f = q + N;

    di = f + N;
    L = di + N;
    U = L + ig[N];

    x = U + ig[N];

    gridX = x + N;
    gridY = gridX + nX;
    edges = gridY + N; 
}