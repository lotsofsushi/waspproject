#include <iostream>
#include "array3.h"

using namespace std;

//These are special containers, that make passing multidimensional arrays to functions easier

double& Array2::operator() (size_t x, size_t y) {
    return data.at(x + y * width);
}

double& Array3::operator() (size_t x, size_t y, size_t z) {
    return data.at(x + y * width + z * width * height);
}

double& Array4::operator() (size_t x, size_t y, size_t z, size_t d) {
    return data.at(x + y * width + z * width * height + d * width * height * depth);
}

void printArray4(const Data& param,Array4& arr, int k, int l)
{
    for (int i = 0; i < param.S+1; i++) {
    for (int j = 0; j < param.P+1; j++) {
    if (j == param.P) cout << arr(i,j,k,l) << endl;
    else cout << arr(i,j,k,l) << " ";
    }}
        cout << endl;
}

void printArray3(const Data& param,Array3& arr, int time)
{
    for (int i = 0; i < param.S+1; i++) {
    for (int j = 0; j < param.P+1; j++) {
    if (j == param.P) cout << arr(i,j,time) << endl;
    else cout << arr(i,j,time) << " ";
    }}
        cout << endl;
}

void printArray2(const Data& param,Array2& arr)
{
    for (int i = 0; i < param.S+1; i++) {
    for (int j = 0; j < param.P+1; j++) {
    if (j == param.P) cout << arr(i,j) << endl;
    else cout << arr(i,j) << " ";
    }}
        cout << endl;
}

void printArray2Special(Array2& arr, int D1, int D2)
{
    for (int i = 0; i < D1+1; i++) {
    for (int j = 0; j < D2+1; j++) {
    if (j == D2) cout << arr(i,j) << endl;
    else cout << arr(i,j) << " ";
    }}
        cout << endl;
}
