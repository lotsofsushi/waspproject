#ifndef BACKWARD_H_INCLUDED
#define BACKWARD_H_INCLUDED
#include "array3.h"
struct Data;

void backwardSimulation (const Data& param, Array3& pi, Array3& newpi, Array2& zi, Array2& newzi, Array3& g1, Array3& newg1, Array3& g2, Array3& newg2, Array3& g3, Array3& newg3, Array3& g4, Array3& newg4, Array3& x, Array3& h, double* y,Array3& mb,Array3& mh,Array3& rho, double& expected, double& expected2, int sim);

void calcFloaters(const Data& param, double* vf, Array3& nest, Array3& pi, Array2& zi, Array3& g2, Array3& vb, Array3& vh, double* floaters,Array3& mb,Array3& mh,Array3& rho, int t);
void calcBreeders(const Data& param, double* vf, Array2& zi,Array3& g1,Array3& g2,Array3& g4, Array3& h, Array3& vb,Array3& mb,Array3& mh,Array3& rho, int t);
void calcHelpers(const Data& param, Array3& h,Array3& g1, Array3& g2, Array3& g3, Array3& g4, double* vf, Array3& vh, Array3& vb,Array3& mb,Array3& mh,Array3& rho, int t);
void calcHelpersFoc(const Data& param, Array3& h,Array3& g2,int choice,double* vf, Array3& vhi, Array3& vh, Array3& vb,Array3& mb,Array3& mh,Array3& rho, int t);


void calcColonies(const Data& param, Array3& h, Array2& zi, Array3& g1,Array3& g2,Array3& g4, Array3& vc,Array3& mb,Array3& mh,Array3& rho,int t);

void calcColDecision(const Data& param, Array3& h,  Array3& g1, int choice, Array3& g2, Array3& g4, Array3& vci, Array3& vc, Array3& mb, Array3& mh, Array3& rho,int t);

void decisionsFloaters(const Data& param,  double* vf, double* Vf, Array3& pi, Array3& newpi, Array2& zi, Array3& g2, Array3& vb, Array3& Vb, Array3& vh, Array3& Vh,Array3& mb,Array3& mh,Array3& rho, int sim);
void decisionsBreeders(const Data& param, double* Vf, Array3& pi, Array3& newpi, Array2& zi, Array3& g2, Array3& vb, Array3& Vb, Array3& vh, Array3& Vh,Array3& mb,Array3& mh,Array3& rho, int sim);
void decisionsHelpers(const Data& param, Array3& h, double* vf, Array3& nest, Array3& pi,  Array3& g1, Array3& g2, Array3& g3, Array3& g4, Array3& newg1, Array3& newg2, Array3& newg3, Array3& newg4, Array3& vb, Array3& vh, Array3& vc,Array3& mb,Array3& mh,Array3& rho, int sim);

//void vFloaters(const Data& param, Array3& finding, Array3& strategy, Array3& vb, Array3& vh, Array3& vf); //calculates rep. values for floaters

//double binomialProb(int n, int k, double p);
//int binomialCoeff(int n, int k);


#endif // BACKWARD_H_INCLUDED
