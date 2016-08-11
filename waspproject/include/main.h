#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED


struct Data
{
 int T;
 int S;
 int P;
 int N;
 int cN;
 double X;

 double y0; //initial density of floaters
 double pF; //probability of finding a nest-site

 double mf; //probability that floater will die during one time step
 double ms0;
 double ms1;
 double mbi;
 double mh1;
 double mh2;
 double mh3;
 double mh4;

 double bf;
 double bp;
 double f0;
 double p0;
 double k;
 double km;

 double r; //relatedness
 double delta; //error level



double convPrec;
bool changeError;
int convAveraging;

bool agg;



};





#endif // MAIN_H_INCLUDED
