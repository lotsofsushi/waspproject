#include <iostream>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <fstream>

#include "forward.h"
#include "main.h"
#include "array3.h"

/*
Forward computation consists of the following steps:
1) Reproduction
2) Building structure
2) Floater decisions
3) Mortality
4) Subordinate staying/leaving decisions
*/


using namespace std;

//Forward simulation
void forwardSimulation(const Data& param, Array3& pi, Array2& zi, Array3& g1, Array3& g2, Array3& g3, Array3& g4, Array3& x, Array3& h,Array3& oldh, double* y,Array3& mb,Array3& mh,Array3& rho, double& produced)
{
    //Clear the x cont. (density of nests)
    for (int t = 0; t < param.T-1; t++)
    {
        for (int s = 0; s < param.S+1; s++)
        {
            for (int p = 0; p < param.P+1; p++)
                x(s,p,t) = 0;
        }
    }

    //Initalize nest and floater densities
    x(0,0,0) = param.X; //Initially all nests are empty with no structure
    y[0] = param.y0;

    //Temporary containers tempx, tempy, xj (hence no need for time dimension)
    Array2 tempx(param.S+1,param.P+1); //Temporary container: density of nests of different size and phase (since its a temporary container it does not need a time dimension)

    for (int s = 0; s < param.S+1; s++)
    {
        for (int p = 0; p < param.P+1; p++) tempx(s,p)=x(s,p,0);
    }


    double tempy = param.y0; //Temporary container for floaters

    Array2 xj(param.S+1,param.P+1); //additive joining effect
    Array4 H(param.S+1,param.P+1,param.S+1,param.P+1); //Building structure probability transition tensor


    produced=0; //set number of produced osspring to 0


    for(int t = 0; t < param.T-1; t++) //t=0 is initialized, so we will start with t=1
    {


        //cout << "tempy: " << tempy << endl;

        //cout << "tempx" << endl;
        //printArray2(param,tempx);


        //cout << " Reproduction!" << endl;
        reproduction(param,pi,zi,g1,g4,x,produced,t);

        //cout << "tempy: " << tempy << endl;

        //cout << " Transformation of nests!" << endl;
        trans(param,zi,g1,g2,g3,g4,x,tempx,tempy,H,xj,rho,mh,mb,t);


        //cout << " Joining nests!" << endl;
        joining(param,y,pi,x,h,oldh,tempy,H,xj,t);


        //update x=L(M(B(x)))+B(xj)
        for(int s = 0; s < param.S+1; s++)
        {
            for(int p = 0; p < param.P+1; p++)
            {
                tempx(s,p)=tempx(s,p)+xj(s,p);
                x(s,p,t+1)=tempx(s,p);
            }
        }

        y[t+1] = tempy;

        //cout << "tempx" << " t " << t <<endl;
        //printArray2(param,tempx);

        //cout << "tempy: " << tempy << endl;

        check(param,y, x,t);


    }

    reproduction(param,pi,zi,g1,g4,x,produced,param.T-1); //reproduction during the the last time step (T-1).
    //cout << "x" << endl;
    //printArray3(param,x,param.T-1);

}

void reproduction(const Data& param, Array3& pi,Array2& zi, Array3& g1, Array3& g4, Array3& x, double& produced, int t)
{
    //Reproduction only happens in nests in the last phase and depends on  the size of the nest and behavior of helpers
    double f[param.S+1]; //reproduction function

    f[0] = 0; //no reproduction in empty nests
    f[1] = param.f0;
    if (param.S > 1)
    {
        for(int s = 2; s < param.S + 1; s++)
        {
            if (param.agg==0) f[s] = param.f0*((1.0-param.k*g1(s,param.P,t))+param.bf*g4(s,param.P,t)*(s-1));
            else f[s] = param.f0*((1.0-param.k*g1(s,param.P,t)*double(s-1)/double(param.S-1))+param.bf*g4(s,param.P,t)*(s-1));

        }

    }


    double  offspring = 0;
    offspring = x(1,param.P,t)*f[1]*(1.0-zi(param.P,t));
    for(int s = 2; s < param.S + 1; s++)
    {
        offspring = offspring + x(s,param.P,t)*f[s];
    }

    produced = produced + offspring;

}

void trans(const Data& param, Array2& zi,Array3& g1, Array3& g2, Array3& g3, Array3& g4, Array3& x, Array2& tempx, double& tempy, Array4& H, Array2& xj,Array3& rho, Array3& mh, Array3& mb, int t)
{

    for (int s=0; s < param.S+1; s++)
    {
        for (int p=0; p < param.P+1; p++)
        {
            for (int news=0; news < param.S+1; news++)
            {
                for (int newp=0; newp < param.P+1; newp++)
                {
                    H(s,p,news,newp)=0;

                }
            }
        }
    }

    Array2 phi(param.S+1,param.P+1);

    //Progression level in each nest
    for (int s = 1; s<param.S+1; s++)
    {
        for (int p = 0; p<param.P; p++)
        {
            if (s == 1) rho(s,p,t) = param.p0;
            else
            {
                if (param.agg==0) rho(s,p,t) = param.p0*((1-param.k*g1(s,p,t)) + param.bp*g4(s,p,t)*(s-1));
                else rho(s,p,t) = param.p0*((1-param.k*g1(s,param.P,t)*double(s-1)/double(param.S-1)) + param.bp*g4(s,p,t)*(s-1));
            }
        }
    }

    //phi
    for (int s = 0; s<param.S+1; s++)
    {
        for (int p = 0; p<param.P+1; p++) phi(s,p)= rho(s,p,t)-floor(rho(s,p,t));
    }

    //mortality
    for(int s = 1; s < param.S+1; s++)
    {
        for(int p = 0; p < param.P+1; p++)
        {
            mh(s,p,t) = g1(s,p,t)*param.mh1 + g2(s,p,t)*param.mh2 + g3(s,p,t)*param.mh3 + g4(s,p,t)*param.mh4;
            if (s==1) mb(s,p,t)= zi(p,t)*param.ms0 + (1-zi(p,t))*param.ms1;
            else
            {
                if (param.agg==0) mb(s,p,t)= param.mbi + param.km*g1(s,p,t);
                else mb(s,p,t)= param.mbi + param.km*g1(s,p,t)*double(s-1)/double(param.S-1);
            }

        }
    }

    //cout << "rho(t=" << t <<")"<< endl;
    //printArray3(param,rho,t);


    //Probability transition tensor s=0
    for (int p=0; p < param.P+1; p++)
    {
        for (int news=0; news < param.S+1; news++)
        {
            for (int newp=0; newp < param.P+1; newp++)
            {
                if (news==0 && p==newp) H(0,p,news,newp) = 1; //transition to the same state
                else   H(0,p,news,newp) = 0;
            }
        }
    }

    //Probability transition tensor s==1
    for (int p=0; p < param.P+1; p++)
    {
        for (int news=0; news < param.S+1; news++)
        {
            for (int newp=0; newp < param.P+1; newp++)
            {
                if (news==1)   //stays and survives
                {

                    if (newp==param.P)   //next phase is last phase
                    {
                        if (floor(rho(1,p,t))+p>=param.P) H(1,p,news,newp)= (1-zi(p,t))*(1-param.ms1); //arrives definitely to the last phase
                        else if (ceil(rho(1,p,t))+p==param.P && rho(1,p,t)!=floor(rho(1,p,t))) H(1,p,news,newp) = phi(1,p)*(1-zi(p,t))*(1-param.ms1); //arrives with phi
                        else H(1,p,news,newp) = 0; //does not arrive to last phase from current phase
                    }

                    else if (newp>=p && newp!=param.P)   //next phase is not last phase
                    {
                        if (ceil(rho(1,p,t))+p==newp && rho(1,p,t)!=floor(rho(1,p,t))) H(1,p,news,newp) = phi(1,p)*(1-zi(p,t))*(1-param.ms1); //upper level progression
                        else if (floor(rho(1,p,t))+p==newp) H(1,p,news,newp) = (1-phi(1,p))*(1-zi(p,t))*(1-param.ms1); //lower level progression
                        else H(1,p,news,newp) = 0;

                    }

                    else H(1,p,news,newp) = 0; //nest phase can not regress

                }

                else if (news==0)   //dies or leaves
                {
                    if (newp==param.P)    //next phase is last phase
                    {
                        if (newp>p)  //progression has taken place
                        {
                            if (floor(rho(1,p,t))+p>=param.P && p!=param.P) H(1,p,news,newp)= (1-zi(p,t))*param.ms1; //arrives definitely to last phase
                            else if (ceil(rho(1,p,t))+p==param.P && rho(1,p,t)!=floor(rho(1,p,t))) H(1,p,news,newp) = phi(1,p)*(1-zi(p,t))*param.ms1; //arrives with phi
                            else H(1,p,news,newp) = 0; //does not arrive to last phase from current phase
                        }

                        else if (newp==p)  H(1,p,news,newp) = (1-zi(p,t))*param.ms1 + zi(p,t); //either stays but dies or leaves
                        else H(1,p,news,newp) = 0;
                    }

                    else if (newp>=p && newp!=param.P)    //next phase is not last phase
                    {
                        if (newp>p)  //progression happens
                        {
                            if (ceil(rho(1,p,t))+p==newp && rho(1,p,t)!=floor(rho(1,p,t))) H(1,p,news,newp) = phi(1,p)*(1-zi(p,t))*param.ms1; //upper level progression
                            else if (floor(rho(1,p,t))+p==newp && floor(rho(1,p,t))>0) H(1,p,news,newp) = (1-phi(1,p))*(1-zi(p,t))*param.ms1; //lower level progression
                            else H(1,p,news,newp) = 0;
                        }
                        else   //progression does not happen
                        {
                            if (floor(rho(1,p,t))==0) H(1,p,news,newp) = (1-phi(1,p))*(1-zi(p,t))*param.ms1+zi(p,t); //tries to progress but fails or leaves
                            else if (floor(rho(1,p,t))!=0) H(1,p,news,newp) = zi(p,t); //leaves
                            else H(1,p,news,newp) = 0;
                        }
                    }

                    else H(1,p,news,newp) = 0; //nest phase can not regress
                }

                else H(1,p,news,newp) = 0; //nest size can not grow
            }
        }
    }



    //Probability transition tensor:multifoundress nests
    for (int s=2; s < param.S+1; s++)
    {
        for (int p=0; p < param.P+1; p++)
        {
            for (int news=0; news < param.S+1; news++)
            {
                for (int newp=0; newp < param.P+1; newp++)
                {

                    //nobody leaves, nobody dies
                    if (news==s)
                    {

                        double prob0 = (1-mb(s,p,t))*bin(s-1,0,mh(s,p,t))*bin(s-1,0,g2(s,p,t)); //probability of nobody dying and leaving
                        //if (newp==0 && p==0) cout << "prob0=" << prob0 << " s=" << s  <<  " news=" << news << endl;
                        if (newp==param.P)
                        {
                            if (floor(rho(s,p,t))+p>=param.P) H(s,p,news,newp)= prob0; //arrives definitely to last phase
                            else if (p==param.P) H(s,p,news,newp)= prob0;
                            else if (ceil(rho(s,p,t))+p==param.P && rho(s,p,t)!=floor(rho(s,p,t))) H(s,p,news,newp) = phi(s,p)*prob0; //arrives with phi
                            else H(s,p,news,newp) = 0; //does not arrive to last phase from current phase
                        }


                        else if (newp>=p && newp!=param.P)   //next phase is not last phase
                        {
                            if (ceil(rho(s,p,t))+p==newp && rho(s,p,t)!=floor(rho(s,p,t))) H(s,p,news,newp) = phi(s,p)*prob0; //upper level progression
                            else if (floor(rho(s,p,t))+p==newp) H(s,p,news,newp) = (1-phi(s,p))*prob0; //lower level progression
                            else H(s,p,news,newp) = 0;
                        }

                        else H(s,p,news,newp) = 0; //nest phase can not regress


                    }

                    //nest size decreases
                    else if (news<s)
                    {
                        double prob1=0;
                        double prob2=0;
                        for (int i=0; i<s-news; i++)
                            prob1 = prob1 + mb(s,p,t)*bin(s-1,i,mh(s,p,t))*bin(s-1-i,news,(1-g2(s,p,t))); //dominant dies
                        for (int i=0; i<s-news+1; i++)
                            prob2 = prob2 + ( news!=0? (1-mb(s,p,t))*bin(s-1,i,mh(s,p,t))*bin(s-1-i,news-1,(1-g2(s,p,t))) :0 ); //dominant survives
                        //if (newp==0 && p==0) cout << "prob1=" << prob1 << " prob2=" << prob2 << " s=" << s  <<  " news=" << news << endl;

                        if (newp==param.P)    //next phase is last phase
                        {
                            if (newp>p)  //progression has taken place
                            {
                                if (floor(rho(s,p,t))+p>=param.P && p!=param.P) H(s,p,news,newp)= prob1+prob2; //arrives definitely to last phase
                                else if (ceil(rho(s,p,t))+p==param.P && rho(s,p,t)!=floor(rho(s,p,t))) H(s,p,news,newp) = phi(s,p)*(prob1+prob2); //arrives with phi
                                else H(s,p,news,newp) = 0; //does not arrive to last phase from current phase
                            }

                            else if (newp==p)  H(s,p,news,newp) = prob1+prob2; //either stays but dies or leaves
                            else H(s,p,news,newp) = 0;
                        }

                        else if (newp>p && newp!=param.P)  //progression happens and not to last phase
                        {
                            if (ceil(rho(s,p,t))+p==newp && rho(s,p,t)!=floor(rho(s,p,t))) H(s,p,news,newp) = phi(s,p)*(prob1+prob2); //upper level progression
                            else if (floor(rho(s,p,t))+p==newp && floor(rho(s,p,t))>0) H(s,p,news,newp) = (1-phi(s,p))*(prob1+prob2); //lower level progression
                            else H(s,p,news,newp) = 0;
                        }

                        else if (newp==p && newp!=param.P)   //progression does not happen
                        {
                            if (floor(rho(s,p,t))==0) H(s,p,news,newp) = (1-phi(s,p))*(prob1+prob2); //tries to progress but fails or leaves
                            else H(s,p,news,newp) = 0;
                        }

                        else H(s,p,news,newp) = 0; //nest phase can not regress


                    }

                    //nest size can not grow and they can not end up s=0
                    else H(s,p,news,newp) = 0;






                }
            }
        }
    }


    /*      for (int s=0; s<param.S+1; s++)
          {
              for (int p=0; p<param.P+1; p++)
              {
                  cout <<"t="<< t << " H(news=" << s << ",newp=" << p <<")"<< endl;
                  printArray4(param,H,s,p);

              }
          }*/


    double leavers=0;

    //solitary breeders who leave
    for(int p = 0; p < param.P+1; p++)
        leavers=leavers+x(1,p,t)*zi(p,t)*(1-param.ms0);

    //subordinates who leave
    for(int s = 2; s < param.S+1; s++)
    {
        for(int p = 0; p < param.P+1; p++)
        {

            double sumleavers = 0;
            for(int i = 0; i < s; i++)
            {
                for(int j = 0; j < s-i; j++)
                {
                    sumleavers = sumleavers+j*bin(s-1,i,mh(s,p,t))*bin(s-1-i,j,g2(s,p,t))*tempx(s,p);
                }
            }
            leavers = leavers + sumleavers;
        }
    }

    //update floaters
    tempy = tempy*(1-param.mf) + leavers;


    Array2 newx(param.S+1,param.P+1);

    //Apply H on x
    for(int news = 0; news < param.S+1; news++)
    {
        for(int newp = 0; newp < param.P+1; newp++)
        {
            for(int s = 0; s < param.S+1; s++)
            {
                for(int p = 0; p < param.P+1; p++) //sum over initial nest phases
                    newx(news,newp) = newx(news,newp) + H(s,p,news,newp)*tempx(s,p);
            }
        }
    }



    //update x with structure building
    for(int s = 0; s < param.S+1; s++)
    {
        for(int p = 0; p < param.P+1; p++)
            tempx(s,p)=newx(s,p);
    }


    //cout << "Hx(t=" << t <<")"<< endl;
    //printArray2(param,tempx);

    //cout << "t=" << t << endl;
    //probSums(param,B,param.P,param.P,param.S,t,0);



    //cout << "leavers=" << leavers << endl;


    /*
               for (int s=1; s < 2; s++)
               {
                   for (int p=0; p < param.P+1; p++)
                   {
                       for (int news=0; news < param.S+1; news++)
                       {
                           for (int newp=0; newp < param.P+1; newp++)
                           {

                cout << "H(" << s << "," << p << "," << news << "," << newp << ")="<< H(s,p,news,newp) << endl;
                }}}}
    */


    //cout << "s= " << s << " p= " << p << " news= " << news << " newp= " << newp << endl;





    /*       for (int s=0; s<param.S+1; s++)
           {
               for (int p=0; p<param.P+1; p++)
               {

                   double sum=0;

                   for (int news=0; news<param.S+1; news++)
                   {
                       for (int newp=0; newp<param.P+1; newp++)
                       {
                           sum=sum+H(s,p,news,newp);
                       }
                   }
                   cout << "t="<< t << " sum(s=" << s << ",p=" << p <<")="<< sum << endl;
               }
           }

    */



}

void joining(const Data& param, double* y, Array3& pi, Array3& x, Array3& h, Array3& oldh, double& tempy, Array4& H,Array2& xj, int t)
{
    //This function calculates the inflow and outflow from nest(s,p,t) by the process of floater joining
    //we add xj(s,p,t) later to tempx(s,p,t)

    Array2 xa(param.S+1,param.P+1); //density of individuals joining nest of type (s,p)
    Array2 newxa(param.S+1,param.P+1);
    Array2 xin(param.S+1,param.P+1); //for storing the inflow into nest of particular type
    Array2 xout(param.S+1,param.P+1); //for storing the inflow into nest of particular type

    double mu = y[t]*param.pF/param.X; //mean number of floaters that are first in current time step to find a nest site
    double joiners = 0; //counting how many join all together


    for (int s=0; s<param.S+1; s++)
    {
        for (int p=0; p<param.P+1; p++)
        {
            oldh(s,p,t) = y[t]*(1.0-param.mf)*param.pF*pi(s,p,t)/param.X;
            h(s,p,t) = (1.0-param.mf)*(1-exp(-1.0*mu))*pi(s,p,t);
            if (h(s,p,t)>1.0) cout << "Error: probability of joining a nest cannot be larger than 1!"  << " Time: " << t << endl;
            if (h(s,p,t)<0.0) cout << "Error: probability of joining a nest cannot be less than 0!"  << endl;
            joiners = joiners+h(s,p,t)*x(s,p,t);
            xa(s,p)=x(s,p,t)*h(s,p,t);
            //cout << "xa(" << s << "," << "p)=" << xa(s,p) << " h= " <<  h(s,p,t) << " x= " << x(s,p,t) << endl;
            //cout << "h[" << s << "]["<< p << "]" << "]["<< t << "]"<< h(s,p,t) << endl;
            //cout << "joiners: " << joiners << " s= " << s << " p= " << p << endl;
            //cout << "oldh[" << s << "]["<< p << "]" << "]["<< t << "]"<< oldh(s,p,t) << endl;
            xj(s,p)=0;
        }
    }


    if (tempy<joiners) cout << "Error: more joiners than floaters!"  << endl;


    //also update how many floaters are left
    tempy = tempy - joiners;


    //Apply H on xin and xout
    for(int news = 0; news < param.S+1; news++)
    {
        for(int newp = 0; newp < param.P+1; newp++)
        {
            for(int s = 0; s < param.S+1; s++)
            {
                for(int p = 0; p < param.P+1; p++)
                {
                    newxa(news,newp) = newxa(news,newp) + H(s,p,news,newp)*xa(s,p);
                }
            }
        }
    }


    for (int p=0; p<param.P+1; p++) //note that floaters who join nest(s,p) end up in nest(s,min(s+p,param.nestPhase)) because of the simultaneous structure building
    {
        xout(0,p) = -newxa(0,p); //nests of size 0 (only outflux)
        xin(param.S,p) = newxa(param.S-1,p); //those who join nests of size S-1 (only influx)
        for (int s=1; s<param.S; s++)  //intermediate sized nests
        {
            xin(s,p) = newxa(s-1,p);
            xout(s,p) = -newxa(s,p);

        }
    }

    //cout << "xa:" << endl;
    //printArray2(param,xa);

    //cout << "newxa:" << endl;
    //printArray2(param,newxa);


    //cout << "xin:" << endl;
    //printArray2(param,xin);

    //cout << "xout:" << endl;
    //printArray2(param,xout);


    for(int s = 0; s < param.S+1; s++)
    {
        for(int p = 0; p < param.P+1; p++)  xj(s,p)=xin(s,p)+xout(s,p);

    }

    //checking if all elements of xj add up to zero
    double sum = 0;
    for (int s=0; s<param.S+1; s++)
    {
        for (int p=0; p<param.P+1; p++) sum = sum + xj(s,p);
    }

    //joiningEffect matrix elements have to add up to 0 or give error message
    if (!almostEqual(sum, 0.0, 0.0000001)) cout << "Error: something wrong with joining nests effect! "  << sum << endl;


}



void check(const Data& param, double* y, Array3& x, int t)
{
    //Check if total nest sites remains constant in time, it should be equal to param.siteDensity across time
    double totalNestSites = 0;
    for (int s = 0; s < param.S+1; s++)
    {
        for (int p = 0; p < param.P+1; p++)
        {
            totalNestSites = totalNestSites + x(s,p,t);
        }
    }
    if (!almostEqual(totalNestSites, param.X, 0.000000001)) cout << "Error: Total nest sites is not constant! "  << endl;

    //Check if the densities are not negative
    if (y[t+1]<0)
        cout << "Error: Negative floater densities at t=" << t+1 << endl;

    for (int s = 0; s < param.S+1; s++)
    {
        for (int p = 0; p < param.P+1; p++)
        {
            if (x(s,p,t)<0)
                cout << "Error: Negative nest densities at t=" << t+1 << endl;

        }

    }



}

void probSums(const Data& param, Array3& arr, int D1, int D2, int D3, int t, int type)
{
    double sum;

    for(int x = 0; x < D1+1; x++)
    {
        for(int z = 0; z < D3+1; z++)
        {
            sum = 0;
            for(int y = 0; y < D2+1; y++)
            {
                sum = sum+arr(x,y,z);
            }
            if (!almostEqual(sum, 1, 0.000000001))
            {
                cout << "Error: Probability transition function ";
                if (type==0) cout << "B" << endl;
                else if (type==1) cout << "M" << endl;
                else if (type==2) cout << "L" << endl;
                else  cout << endl;
                cout << "sum = " << sum << "("<< x << "," << z << ") t=" << t << endl;

                printArray3(param,arr,z);

            }
        }
    }

}

double bin(int n, int k, double p)
{
    if (k>n) cout << "Error! Binomial: k>n!" << endl;
    double res;
    res = binomialCoeff(n,k)*pow(p,k)*pow(1-p,n-k);
    return res;

}

int binomialCoeff(int n, int k)
{
    int res = 1;

    // Since C(n, k) = C(n, n-k)
    if ( k > n - k ) k = n - k;

    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }

    return res;
}

bool almostEqual(double a, double b, double delta)
{
    double const diff = std::abs (a - b);
    double const toll = delta;
    //std::numeric_limits<double>::epsilon();
    return diff < toll;
}







