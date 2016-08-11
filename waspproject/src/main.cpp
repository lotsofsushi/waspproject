#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>

#include "main.h"
#include "array3.h"
#include "forward.h"
#include "backward.h"
#include "print.h"
#include "init.h"

/*
03/08 version

Dynamic programming model for describing semisociality in wasps.

State variables: s and p
Floaters and subordinates can make decisions
Limitations of the implementation:
Nest size can not be larger than 25 individuals, since the binomial function is implemented in such a way that
int cannot hold that high values.

*/

using namespace std;

//void mortality (double arr[][nestSize+1]);
//void buildingStructure (double arr[][nestSize+1]);
//void joiningNests (double arr[][nestSize+1], double floaters, double y, double pi[][nestSize+1],double nestSites,int t);


int main()
{
    //Parameter values:
    Data param;
    param.T = 70; //number of time steps
    param.S = 6; //max nest size
    param.P = 4; //max nest phase
    param.N = 400;
    param.cN = 200;

    param.X = 400; //max density of possible nests site
    param.y0 = 400; //density of floaters
    param.pF = 0.8; //probability of finding a nest-site

    param.mf = 0.008; //probability that floater will die during one time step
    param.ms0 = 0.05; //probability to die of a solitary female if she leaves
    param.ms1 = 0.05; //probability to die of a solitary female if she stays
    param.mbi = 0; //probability that dominant female, that has at least one subordinate
    param.mh1 = 0.008; //probability of death, being agressive subordinate
    param.mh2 = 0.008; //probability of death, leaving the nest
    param.mh3 = 0.008; //probability of death, sitting in the nest
    param.mh4 = 0.015; //probability of death, foraging

    param.f0 = 10; //intrinsic fecundity
    param.p0 = 0.1; // intrinsic progression rate
    param.bf = 1; //fecundity increase by a helper
    param.bp = 1; //colony progression increase

    param.convPrec = 0.05; //Convergence precision
    param.changeError = 1; //0-don't wanna see, 1- show
    param.convAveraging = 4; //how many iterations of averaging

    param.k = 0.5; // breeders susceptibility to aggression (effect on productivity) 0...1
    param.km = 0.01;//1-param.mbi; // breeders susceptibility to aggression (effect on mortality) 0...(1-mbi)
    param.agg = 1; //0 - old way, 1- new way


    param.r = 0.25; //relatedness
    param.delta = 2; //delta>0, controls for the degree of error, if delta->0 then the probability of making mistakes tend to 0

    double expected; //expected offspring
    double expected2; //expected offspring
    double produced; //produced offspring

    Array2 change(param.N,5); //sum of all differences of elements in a decision matrix


    Array3 x(param.S+1,param.P+1,param.T); //density of nests of different size and phase at each time step

    //Decision matrices
    Array3 pi(param.S+1,param.P+1,param.T); //prob. to join the nest (s,p) at t
    Array2 zi(param.P+1,param.T); //prob. to leave the nest, of a solitary female
    Array3 g1(param.S+1,param.P+1,param.T); //prob. to be aggressive in the nest (s,p) at t
    Array3 g2(param.S+1,param.P+1,param.T); //prob. to leave the nest (s,p) at t
    Array3 g3(param.S+1,param.P+1,param.T); //prob. to sit in the nest (s,p) at t
    Array3 g4(param.S+1,param.P+1,param.T); //prob. to forage in the nest (s,p) at t

    //New decision matrices
    Array3 newpi(param.S+1,param.P+1,param.T);
    Array2 newzi(param.P+1,param.T);
    Array3 newg1(param.S+1,param.P+1,param.T);
    Array3 newg2(param.S+1,param.P+1,param.T);
    Array3 newg3(param.S+1,param.P+1,param.T);
    Array3 newg4(param.S+1,param.P+1,param.T);

    Array3 h(param.S+1,param.P+1,param.T); //probability that the nest will be joined by a helper at time t
    Array3 oldh(param.S+1,param.P+1,param.T); //probability that the nest will be joined by a helper at time t

    //Mortality
    Array3 mb(param.S+1,param.P+1,param.T);
    Array3 mh(param.S+1,param.P+1,param.T);

    //Building structure function
    Array3 rho(param.S+1,param.P+1,param.T);

    //Floaters
    double floaters[param.T];
    double* y = &floaters[0];

    printParam(param);

    //Initial conditions for strategies
    init_pi(param,pi, 0.9, 0.1); //floater decisions
    init_zi(param,zi,0);//leaving the solitary nest
    init_gi(param,g1,g2,g3,g4,0.25,0.25,0.25,0.25); //helper decisions

    for (int n = 0; n < param.N; n++) {

    forwardSimulation(param,pi,zi,g1,g2,g3,g4,x,h, oldh,y,mb,mh,rho, produced);
    backwardSimulation(param,pi,newpi,zi,newzi,g1,newg1,g2,newg2,g3,newg3,g4,newg4,x,h,y,mb,mh,rho, expected,expected2,n);


    for(int t=0;t<param.T;t++){
    for(int s=0;s<param.S+1;s++){
    for(int p=0;p<param.P+1;p++){
    change(n,0) = change(n,0) + abs(pi(s,p,t)-newpi(s,p,t));
    change(n,1) = change(n,1) + abs(g1(s,p,t)-newg1(s,p,t));
    change(n,2) = change(n,2) + abs(g2(s,p,t)-newg2(s,p,t));
    change(n,3) = change(n,3) + abs(g3(s,p,t)-newg3(s,p,t));
    change(n,4) = change(n,4) + abs(g4(s,p,t)-newg4(s,p,t));


    if (n==param.N-1 && param.changeError==1) {
        if(abs(pi(s,p,t)-newpi(s,p,t))>param.convPrec)
        cout << "pi(" << s << "," << p << "," << t << "). Change="<<  abs(pi(s,p,t)-newpi(s,p,t)) << endl;
        if(abs(g1(s,p,t)-newg1(s,p,t))>param.convPrec)
        cout << "g1(" << s << "," << p << "," << t << "). Change="<<  abs(g1(s,p,t)-newg1(s,p,t)) << endl;
        if(abs(g2(s,p,t)-newg2(s,p,t))>param.convPrec)
        cout << "g2(" << s << "," << p << "," << t << "). Change="<<  abs(g2(s,p,t)-newg2(s,p,t)) << endl;
        if(abs(g3(s,p,t)-newg3(s,p,t))>param.convPrec)
        cout << "g3(" << s << "," << p << "," << t << "). Change="<<  abs(g3(s,p,t)-newg3(s,p,t)) << endl;
        if(abs(g4(s,p,t)-newg4(s,p,t))>param.convPrec)
        cout << "g4(" << s << "," << p << "," << t << "). Change="<<  abs(g4(s,p,t)-newg4(s,p,t)) << endl;
    }

    pi(s,p,t) = newpi(s,p,t);
    g1(s,p,t) = newg1(s,p,t);
    g2(s,p,t) = newg2(s,p,t);
    g3(s,p,t) = newg3(s,p,t);
    g4(s,p,t) = newg4(s,p,t);


    }}}

    //if(n % 100 == 0) cout << n << endl;
    cout << "Change(" << n << ")="<< change(n,0) << "; " << change(n,1) << "; " << change(n,2) << "; " << change(n,3) << "; " << change(n,4) <<endl;
    cout << "produced[" << n << "]="<< produced << endl;
    cout << "expected[" << n << "]="<< expected << endl;
    cout << "expected2[" << n << "]="<< expected2 << endl;
    cout << "error = "<< abs(produced-expected)/abs(expected)*100 << "%"<< endl;
    cout << "error = "<< abs(produced-expected2)/abs(expected2)*100 << "%"<< endl;
    //for (int t = 0; t < param.T-1; t++) {
    //cout << "pi matrix n = " << n << " and t = " << t << endl;
    //printArray3(param,pi,t);}

    //printArray2Special(zi,param.P,param.T-1);

    if (n==param.N-1)
    printData(param,x,pi,g1,g2,g3,g4,y);

    }




    return 0;
}





