#ifndef COMMON_H
#define COMMON_H

#include <fstream>
#include <iostream>

extern int nX, nY, N; 
extern int numIter, i, k, step;
extern int *ig, *jg;
extern double *x, *f;
extern double *gridX, *gridY;
extern double *edges;
extern double *L, *U, *ggl, *ggu; 
extern double *di, *diOrig;
extern double *r, *p, *z, *q;
extern double eps;

void readSize();
void readGrid();
void config();
void buildPortrait();
void readEdges();
double getAnalytical(double, double);
double getF(double, double);
double getLambda(double, double);
double getGamma(double, double);
double getBeta();
void globalAsm();
void boundaryCond();
void firstCond(int);
void secondCond(int);
void thirdCond(int);
double sum(int i, int j);
void factorization();
void addGlobal(double, int, int);
void forwardGauss(double *in, double *out);
void backwardGauss(double *in, double *out);
void multVectByA(double *in, double *out);
double multVectByVect(double *v1, double *v2);
void los();
void printResult();

#endif
