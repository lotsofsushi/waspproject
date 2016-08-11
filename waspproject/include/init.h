#ifndef INIT_H_INCLUDED
#define INIT_H_INCLUDED

#include <iostream>
#include "main.h"
#include "array3.h"

void init_pi(const Data& param, Array3& pi, double s0, double si);
void init_zi(const Data& param, Array2& zi, double z0);
void init_gi(const Data& param, Array3& g1, Array3& g2, Array3& g3, Array3& g4, double g10, double g20, double g30, double g40);


#endif // INIT_H_INCLUDED
