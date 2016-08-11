#ifndef ARRAY3_H_INCLUDED
#define ARRAY3_H_INCLUDED

#include <iostream>
#include <vector>
#include "main.h"

class Array2 {
    size_t width, height;
    std::vector<double> data;
  public:
    Array2(size_t x, size_t y, int init = 0):
      width(x), height(y), data(x*y, init) {}

    double& operator() (size_t x, size_t y);
};


class Array3 {
    size_t width, height;
    std::vector<double> data;
  public:
    Array3(size_t x, size_t y, size_t z, int init = 0):
      width(x), height(y), data(x*y*z, init) {}

    double& operator() (size_t x, size_t y, size_t z);
};

class Array4 {
    size_t width, height, depth;
    std::vector<double> data;
  public:
    Array4(size_t x, size_t y, size_t z, size_t d, int init = 0):
      width(x), height(y), depth(z), data(x*y*z*d, init) {}

    double& operator() (size_t x, size_t y, size_t z, size_t d);
};
// Kasutus
// Array3 arr(10, 10, 10); // 10x10x10 array, mille iga liige on "0"
// arr(0, 0, 0) = 1; // esimese liikme väärtus on nüüd "1"

void printArray4(const Data& param,Array4& arr, int k, int l);
void printArray3(const Data& param, Array3& arr, int time);
void printArray2(const Data& param, Array2& arr);
void printArray2Special(Array2& arr, int D1, int D2);

#endif // ARRAY3_H_INCLUDED
