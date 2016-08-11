#include <iostream>
#include "init.h"

void init_pi(const Data& param, Array3& pi, double s0, double si) {

    for (int t = 0; t < param.T; t++) {
    for (int p = 0; p < param.P + 1; p++){
    pi(0,p,t) = s0;
    pi(param.S,p,t) = 0; //females will not join nests with maximum size
    for (int s = 1; s < param.S; s++) {
    pi(s,p,t) = si;
    }}}

}

void init_zi(const Data& param, Array2& zi, double z0) {

    for (int t = 0; t < param.T; t++) {
    for (int p = 0; p < param.P + 1; p++){
    zi(p,t) = z0;
    }}

}

void init_gi(const Data& param, Array3& g1, Array3& g2, Array3& g3, Array3& g4, double g10, double g20, double g30, double g40) {


    if (param.S>1)
    {
        for (int s = 2; s < param.S+1; s++) { //helpers are only present in colonies of size 2 and larger
        for (int t = 0; t < param.T; t++) {
        for (int p = 0; p < param.P + 1; p++){
            g1(s,p,t) = g10;
            g2(s,p,t) = g20;
            g3(s,p,t) = g30;
            g4(s,p,t) = g40;
        }}}


    }


}
