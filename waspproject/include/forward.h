#ifndef JOINING_H_INCLUDED
#define JOINING_H_INCLUDED
#include "array3.h"
struct Data;

void forwardSimulation (const Data& param, Array3& pi, Array2& zi, Array3& g1, Array3& g2, Array3& g3, Array3& g4, Array3& nest, Array3& h, Array3& oldh,double* floaters,Array3& mb,Array3& mh,Array3& rho, double& produced);

void reproduction(const Data& param, Array3& pi, Array2& zi, Array3& g1, Array3& g4, Array3& x, double& produced, int t);
void trans(const Data& param, Array2& zi, Array3& g1,Array3& g2, Array3& g3, Array3& g4, Array3& x, Array2& tempx, double& tempy,Array4& H, Array2& xj, Array3& rho, Array3& mh, Array3& mb, int t);
void joining(const Data& param, double* y, Array3& pi, Array3& x, Array3& h, Array3& oldh, double& tempy, Array4& H,Array2& xj, int t);


void check(const Data& param, double* y, Array3& x, int t);
void probSums(const Data& param, Array3& arr, int D1, int D2, int D3, int t, int type);

int binomialCoeff(int n, int k);
double bin(int n, int k, double p);
bool almostEqual ( double a, double b, double delta );




#endif // JOINING_H_INCLUDED
