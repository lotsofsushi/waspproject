#ifndef PRINT_H_INCLUDED
#define PRINT_H_INCLUDED

#include <iostream>
#include "main.h"
#include "array3.h"

void printParam(const Data& param);
void printData(const Data& param, Array3& nest, Array3& pi, Array3& g1, Array3& g2, Array3& g3, Array3& g4, double* floaters);
void printFitness(const Data& param, Array3& vh1, Array3& vh2, Array3& vh3, Array3& vh4, Array3& Vh1, Array3& Vh2, Array3& Vh3, Array3& Vh4, Array3& g1, Array3& g2, Array3& g3, Array3& g4);

#endif // PRINT_H_INCLUDED
