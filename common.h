#ifndef COMMON_H
#define COMMON_H

typedef double real;

extern int nX, nY, N; 
extern int *port, *pos;
extern real *gridX, *gridY;
extern real *L, *U, *Lo, *Uo; 
extern real *diCond, *diOrig;

void readSize();
void readGrid();
void config();
real analytical(real, real);
real getF(real, real);
real getLambda(real, real);
real getGamma(real, real);
//void globalAsm();
//void boundaryCond();
real dotProduct(int i, int j);
int factorization();
void addGlobal(int, int, real);
void lowerGauss(real *in, real *out);
void upperGauss(real *in, real *out);
void multMatVect(real *in, real *out);
real multVect(real *v1, real *v2);
//void result();

#endif
