#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <fstream>


#include "backward.h"
#include "main.h"
#include "array3.h"
#include "forward.h"
#include "print.h"

using namespace std;

/*
Backward computation
Consider a focal individual carrying a rare mutant allele, which is r times as likely
to be present in offspring of other colony members, compared to its own. We
define inclusive reproductive value V(t) as the total number of offspring produced
(at time t, and afterwards) because of the focal individual’s actions, whererelatives’ offspring are discounted by r. We define direct reproductive value v(t)
as the focal individual’s total number of own offspring (produced at time t, and
afterwards). We define colony reproductive value v colony (t) as the total number of
offspring produced by the focal colony (at time t, and afterwards).
Final time, T
The last bout of reproduction happens at T. At this time, direct reproductive
values of y, helpers, and breeders are, respectively,
v f (T) = 0
v h (T) = 0
v b (T,n,w) = f(n,w),
and a colony’s reproductive value is
v colony (T,n,w)= f(n,w).
For y, it is too late at this point to influence anyone’s reproduction, so their
inclusive reproductive value is V f (T) = 0.
For helpers,
V h (T)=[f(n,w)-f(n-1,w)] * r,
which represents the additional offspring produced by the colony owing to the
focal individual’s presence, discounted by relatedness. For breeders,
V b (T)= f(n,w)-f(n-1,w) * r,
which represents the breeder’s direct reproduction, minus the production of
mutant alleles that would have occurred anyway, even in the absence of the
breeder (e.g if the breeder had died, or had failed to join the colony in the
previous time step).
*/


void backwardSimulation (const Data& param, Array3& pi, Array3& newpi,  Array2& zi, Array2& newzi, Array3& g1, Array3& newg1, Array3& g2, Array3& newg2, Array3& g3, Array3& newg3, Array3& g4, Array3& newg4, Array3& x, Array3& h, double* y, Array3& mb,Array3& mh,Array3& rho, double& expected, double& expected2, int sim)
{

    double floaterRV_inclusive[param.T];
    double* Vf = &floaterRV_inclusive[0]; //inclusive RV of y (pointer)
    double floaterRV_direct[param.T];
    double* vf = &floaterRV_direct[0]; //direct RV of y (pointer)

    Array3 Vb(param.S+1,param.P+1,param.T); //inclusive RV of breeders
    Array3 Vh(param.S+1,param.P+1,param.T); //inclusive RV of helpers
    Array3 vb(param.S+1,param.P+1,param.T); //direct RV of breeders (vb(0,p,t) = undefined)
    Array3 vh(param.S+1,param.P+1,param.T); //direct RV of helpers
    Array3 vc(param.S+1,param.P+1,param.T); //direct RV of colonies


    ///////////////////////////////////////////////////////////
    ////////////////Final time/////////////////////////////////
    ///////////////////////////////////////////////////////////

    //RV of y
    vf[param.T-1] = 0;
    Vf[param.T-1] = 0;

    double f[param.S+1]; //reproduction function

    f[0] = 0; //no reproduction in empty nests
    f[1] = param.f0;
    if (param.S > 1)
    {
        for(int s = 2; s < param.S + 1; s++)
        {
            if (param.agg==0) f[s] = param.f0*((1.0-param.k*g1(s,param.P,param.T-1))+param.bf*g4(s,param.P,param.T-1)*(s-1));
            else f[s] = param.f0*((1.0-param.k*g1(s,param.P,param.T-1)*double(s-1)/double(param.S-1))+param.bf*g4(s,param.P,param.T-1)*(s-1));
        }
    }

    //RV of breeders, helpers, colonies

    for(int s = 0; s < 2; s++)
        for(int p = 0; p < param.P + 1; p++)
        {
            vb(s,p,param.T-1) = (s==1?(1-zi(param.P,param.T-1)):1)*(p==param.P?f[s]:0); //if s=0 (there are no breeders in x so 0), if s=c (when p=P)
            Vb(s,p,param.T-1) = vb(s,p,param.T-1);
            vc(s,p,param.T-1) = vb(s,p,param.T-1);
            vh(s,p,param.T-1) = 0; //there are no helpers in colonies of size 0 or 1, so this term really has no meaning
            Vh(s,p,param.T-1) = 0; // -''-
        }

    if (param.S>1)
    {
        for (int s = 2; s < param.S + 1; s++)
        {
            for(int p = 0; p < param.P + 1; p++)
            {
                vb(s,p,param.T-1) = (p==param.P?f[s]:0);
                Vb(s,p,param.T-1) = vb(s,p,param.T-1)-vb(s-1,p,param.T-1)*param.r;
                vc(s,p,param.T-1) = vb(s,p,param.T-1);
                vh(s,p,param.T-1) = 0; //helpers have no direct terminal RV
                Vh(s,p,param.T-1) = (vb(s,p,param.T-1)-vb(s-1,p,param.T-1))*param.r;

            }
        }
    }

    /*
        cout << endl << "vh["<<param.T-1<<"]" << endl;
        printArray3(param,vh,param.T-1);

        cout << endl << "vb["<<param.T-1<<"]" << endl;
        printArray3(param,vb,param.T-1);

        cout << endl << "vc["<<param.T-1<<"]" << endl;
        printArray3(param,vc,param.T-1);

        cout << endl << "Vh["<<param.T-1<<"]" << endl;
        printArray3(param,Vh,param.T-1);

        cout << endl << "Vb["<<param.T-1<<"]" << endl;
        printArray3(param,Vb,param.T-1);


    */

    ///////////////////////////////////////////////////////////
    ////////////////Time t<T///////////////////////////////////
    ///////////////////////////////////////////////////////////
    for (int t=param.T-2; t>=0; t--)
    {


        //y: inclusive reproductive value
        calcFloaters(param,vf,x,pi,zi,g2,vb,vh,y,mb,mh,rho,t);
        calcFloaters(param,Vf,x,pi,zi,g2,Vb,Vh,y,mb,mh,rho,t);

        //Direct RV of Breeders
        calcBreeders(param,vf,zi,g1,g2,g4,h,vb,mb,mh,rho,t);

        //Direct RV of Helpers
        calcHelpers(param,h,g1,g2,g3,g4,vf,vh,vb,mb,mh,rho,t);

        //Direct RV of Colonies
        calcColonies(param,h,zi,g1,g2,g4,vc,mb,mh,rho,t);

        //cout << endl << "vh["<<t<<"]" << endl;
        //printArray3(param,vh,t);
        //cout << endl << "vb["<<t<<"]" << endl;
        //printArray3(param,vb,t);
        //cout << endl << "vc["<<t<<"]" << endl;
        //printArray3(param,vc,t);


        //Inclusive helpers
        for(int s = 2; s < param.S + 1; s++)
        {
            for(int p = 0; p < param.P + 1; p++)
                Vh(s,p,t)=vh(s,p,t) + (vc(s,p,t)-vc(s-1,p,t)-vh(s,p,t))*param.r;
        }

        //cout << endl << "Vh["<<t<<"]" << endl;
        //printArray3(param,Vh,t);

        //cout << endl << "vh["<<t<<"]" << endl;
        //printArray3(param,vh,t);

        //Inclusive breeders
        for(int s = 1; s < param.S + 1; s++)
        {
            for(int p = 0; p < param.P + 1; p++) Vb(s,p,t)=vb(s,p,t) + (vc(s,p,t)-vc(s-1,p,t)-vb(s,p,t))*param.r;
        }

        //cout << endl << "Vb["<<t<<"]" << endl;
        //printArray3(param,Vb,t);


        //cout << endl << "vb["<<t<<"]" << endl;
        //printArray3(param,vb,t);






    } //time loop

    decisionsFloaters(param,vf,Vf,pi,newpi,zi,g2,vb,Vb,vh,Vh,mb,mh,rho,sim);
    decisionsHelpers(param,h,vf,x,pi,g1,g2,g3,g4,newg1,newg2,newg3,newg4,vb,vh,vc,mb,mh,rho,sim);

    /*

    for (int t = 0; t < param.T-1; t++) {
    cout << "v_stay t = " << t << endl;
    printArray3(param,Vs,t);
    cout << "v_leave t = " << t+1 << " is "<< Vl[t] << endl;
    }*/


    expected = Vf[0]*y[0];
    expected2 = vc(0,0,0)*param.X;



}

void calcFloaters(const Data& param, double* vf, Array3& x, Array3& pi, Array2& zi, Array3& g2, Array3& vb, Array3& vh, double* y,Array3& mb,Array3& mh,Array3& rho, int t)
{

    Array2 join(param.S+1,param.P+1); //probability to find and join a nest (s,p)
    Array2 notjoin(param.S+1,param.P+1); //probability to find and not join a nest (s,p)
    Array2 phi(param.S+1,param.P+1); //nest progression probability


    double mu = y[t]*param.pF/param.X; //mean number of y that are first in current time step to find a x site

    for(int s = 0; s < param.S + 1; s++)
    {
        for(int p = 0; p < param.P + 1; p++)
        {
            join(s,p) = (1.0-exp(-1.0*mu))*x(s,p,t)*pi(s,p,t)/y[t]; //find the x site and join
            notjoin(s,p) = (1.0-exp(-1.0*mu))*x(s,p,t)*(1-pi(s,p,t))/y[t]; //find the x site but not join
            phi(s,p)=rho(s,p,t)-floor(rho(s,p,t));
        }
    }

//for(int s = 0; s < param.S + 1; s++){
//for(int p = 0; p < param.P + 1; p++) //cout << "join(" << s << "," << p <<"): "<<  join(s,p) << endl;
//}


    double fTerm[9]; //terms in the floater reproductive value equation
    memset(fTerm, 0, 9*sizeof(double)); //initialize all 0
    double sum1, sum2, sum3;

    for(int p = 0; p < param.P + 1; p++)  //joins nests of size 0 and 1
    {
        fTerm[0] = fTerm[0] + join(0,p)*vb(1,p,t+1); //nests of size 0
        fTerm[1] = fTerm[1] + join(1,p)*(1-zi(p,t))*param.ms1     * (phi(1,p)*vb(1,min(param.P,int(p+ceil(rho(1,p,t)))),t+1)   + (1-phi(1,p))*vb(1,min(param.P,int(p+floor(rho(1,p,t)))),t+1)); //breeders decides to stay but dies
        fTerm[2] = fTerm[2] + join(1,p)*(1-zi(p,t))*(1-param.ms1) * (phi(1,p)*vh(2,min(param.P,int(p+ceil(rho(1,p,t)))),t+1)   + (1-phi(1,p))*vh(2,min(param.P,int(p+floor(rho(1,p,t)))),t+1)); //breeder decides to stay and survives
        fTerm[3] = fTerm[3] + join(1,p)*zi(p,t)*vb(1,p,t+1); //breeder decides to leave

    }

    for(int s = 2; s < param.S; s++)  //goes up to S-1, since floaters can not join x with max size
    {
        for(int p = 0; p < param.P + 1; p++)
        {

            sum1 = 0;
            sum2 = 0;
            sum3 = 0;
            for(int j = 0; j < s; j++)
            {
                for(int r = 0; r < s-j; r++)
                {
                    double ph = bin(s-1,j,mh(s,p,t))*bin(s-1-j,r,g2(s,p,t));
                    sum1 = sum1 + ph * (1.0/double(s-j-r))             * (phi(s,p)*vb(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)   + (1-phi(s,p))*vb(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum2 = sum2 + ph * (double(s-j-r-1)/double(s-j-r)) * (phi(s,p)*vh(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)   + (1-phi(s,p))*vh(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum3 = sum3 + ph *                                   (phi(s,p)*vh(s-j-r+1,min(param.P,int(p+ceil(rho(s,p,t)))),t+1) + (1-phi(s,p))*vh(s-j-r+1,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                }
            }
            fTerm[4] = fTerm[4] + join(s,p)*mb(s,p,t)*sum1;
            fTerm[5] = fTerm[5] + join(s,p)*mb(s,p,t)*sum2;
            fTerm[6] = fTerm[6] + join(s,p)*(1-mb(s,p,t))*sum3;
        }
    }

    for(int s = 0; s < param.S+1; s++)  //goes up to S
    {
        for(int p = 0; p < param.P + 1; p++)
        {
            fTerm[7] = fTerm[7] + notjoin(s,p)*vf[t+1];
        }
    }

    fTerm[8] = (1.0-param.X*(1.0-exp(-1.0*mu))/y[t])*vf[t+1];
    //cout << "Not join=" << (1.0-param.X*(1.0-exp(-1.0*mu))/y[t]) << endl;
    // cout << "Not find=" << (1.0-param.X*(1.0-exp(-1.0*mu))/y[t]) << endl;
    //calculate floater reproductive value for current time step
    //cout << endl << "RV of floaters t= " << t << endl;
    //for (int i=0; i<9; i++) cout << "fterm["<< i<<"]: " << fTerm[i] << " ";
    vf[t] = (1-param.mf)*(fTerm[0]+fTerm[1]+fTerm[2]+fTerm[3]+fTerm[4]+fTerm[5]+fTerm[6]+fTerm[7]+fTerm[8]);
    //cout <<  endl << "vf[" << t << "]: " << vf[t] << endl;
}


void calcBreeders(const Data& param, double* vf,Array2& zi,Array3& g1,Array3& g2,Array3& g4, Array3& h, Array3& vb,Array3& mb,Array3& mh,Array3& rho, int t)
{

    double fTerm[4]; //terms in the floater reproductive value equation
    double sum1, sum2;

    Array2 phi(param.S+1,param.P+1); //nest progression probability
    for(int s = 0; s < param.S + 1; s++)
    {
        for(int p = 0; p < param.P + 1; p++)
            phi(s,p)=rho(s,p,t)-floor(rho(s,p,t));
    }

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

    //Breeders in colonies of size 1
    for(int p = 0; p < param.P + 1; p++)
    {
        fTerm[0] = (1-h(1,p,t))*(phi(1,p)*vb(1,min(param.P,int(p+ceil(rho(1,p,t)))),t+1) + (1-phi(1,p))*vb(1,min(param.P,int(p+floor(rho(1,p,t)))),t+1));
        fTerm[1] = h(1,p,t)*(phi(1,p)*vb(2,min(param.P,int(p+ceil(rho(1,p,t)))),t+1) + (1-phi(1,p))*vb(2,min(param.P,int(p+floor(rho(1,p,t)))),t+1));

        vb(1,p,t) = (1-zi(p,t))*(p==param.P?f[1]:0) + (1-zi(p,t))*(1-param.ms1)*(fTerm[0] + fTerm[1]) + zi(p,t)*(1-param.ms0)*vf[t+1];
        //cout << " t= " << t << " p= " << p << " t0= " << fTerm[0]  << " t1= " << fTerm[1]  << endl;
    }

    //Breeders in colonies of size 2 and larger
    for (int s = 2; s < param.S + 1; s++)
    {
        for(int p = 0; p < param.P + 1; p++)
        {
            sum1=0;
            sum2=0;
            for(int j = 0; j < s; j++)
            {
                for(int r = 0; r < s-j; r++)
                {
                    double ph = bin(s-1,j,mh(s,p,t))*bin(s-1-j,r,g2(s,p,t));
                    sum1 = sum1 + ph * (phi(s,p)*vb(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)                + (1-phi(s,p))*vb(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum2 = sum2 + ph * (phi(s,p)*vb(min(param.S,s-j-r+1),min(param.P,int(p+ceil(rho(s,p,t)))),t+1) + (1-phi(s,p))*vb(min(param.S,s-j-r+1),min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                }
            }
            fTerm[2] = (1-h(s,p,t))*sum1;
            fTerm[3] = h(s,p,t)*sum2;
            //if (p == param.P) cout << fTerm[2] << "|" << fTerm[3] << endl;
            //else cout << fTerm[2] << "|" << fTerm[3] << " ";
            vb(s,p,t) = (p==param.P?f[s]:0) + (1-mb(s,p,t))*(fTerm[2] + fTerm[3]);
        }
    }

}

void calcHelpers(const Data& param, Array3& h,Array3& g1, Array3& g2, Array3& g3, Array3& g4, double* vf, Array3& vh, Array3& vb,Array3& mb,Array3& mh,Array3& rho, int t)
{


    double fTerm[7]; //terms in the helper reproductive value equation
    double act[3]; //probability of choosing actions with appropriate mortality, except leaving the nest
    double sum1, sum2, sum3, sum4, sum5, sum6;

    double g1_; //adjusted value of g1
    double phi_; //adjisted value of phi
    double rho_; //adjusted value of rho
    double mb_; //adjusted value of mb


    for (int s = 2; s < param.S + 1; s++) //there are no helpers in colonies of size 0 and 1 (so they are set to 0 by default)
    {
        for(int p = 0; p < param.P + 1; p++)
        {

            memset(fTerm, 0, 7*sizeof(double)); //set all of them to 0
            act[0] = g1(s,p,t)*(1-param.mh1);
            act[1] = g3(s,p,t)*(1-param.mh3);
            act[2] = g4(s,p,t)*(1-param.mh4);

//
            for (int u = 1; u < 4; u++)   // u-action of the focal individual, u2 is not counted, so sitting is u=2 and foraging is u=3 here
            {

                g1_= (u==1?g1(s,p,t)*double(s-1)/double(param.S-1):g1(s,p,t)*double(s-2)/double(param.S-1));
                if (param.agg==0)
                {
                    rho_ = param.p0*(1-param.k*g1(s,p,t) + param.bp*(g4(s,p,t)*(s-2)+(u==3?1:0)));
                    mb_  = param.mbi + param.km*g1(s,p,t);
                }
                else
                {
                    //rho_ = param.p0*(1-param.k*g1_ + param.b*(g4(s,p,t)*(s-2)+(u==3?1:0)));
                    //mb_= param.mbi + param.km*g1_;
                    rho_=rho(s,p,t);
                    mb_=mb(s,p,t);
                }
                phi_ = rho_-floor(rho_);


                sum1=0;
                sum2=0;
                sum3=0;
                sum4=0;
                sum5=0;
                sum6=0;
                for(int j = 0; j < s-1; j++)
                {
                    for(int r = 0; r < s-j-1; r++)
                    {
                        double ph = bin(s-2,j,mh(s,p,t))*bin(s-2-j,r,g2(s,p,t));
                        sum1 = sum1 + ph * (1.0/double(s-j-r-1))             * (phi_*vb(s-j-r-1,min(param.P,int(p+ceil(rho_))),t+1)              + (1-phi_)*vb(s-j-r-1,min(param.P,int(p+floor(rho_))),t+1));
                        sum2 = sum2 + ph * (1.0/double(s-j-r))               * (phi_*vb(s-j-r,min(param.P,int(p+ceil(rho_))),t+1)                + (1-phi_)*vb(s-j-r,min(param.P,int(p+floor(rho_))),t+1));
                        sum3 = sum3 + ph * (double(s-j-r-2)/double(s-j-r-1)) * (phi_*vh(s-j-r-1,min(param.P,int(p+ceil(rho_))),t+1)              + (1-phi_)*vh(s-j-r-1,min(param.P,int(p+floor(rho_))),t+1));
                        sum4 = sum4 + ph * (double(s-j-r-1)/double(s-j-r))   * (phi_*vh(s-j-r,min(param.P,int(p+ceil(rho_))),t+1)                + (1-phi_)*vh(s-j-r,min(param.P,int(p+floor(rho_))),t+1));
                        sum5 = sum5 + ph *                                     (phi_*vh(s-j-r,min(param.P,int(p+ceil(rho_))),t+1)                + (1-phi_)*vh(s-j-r,min(param.P,int(p+floor(rho_))),t+1));
                        sum6 = sum6 + ph *                                     (phi_*vh(min(param.S,s-j-r+1),min(param.P,int(p+ceil(rho_))),t+1) + (1-phi_)*vh(min(param.S,s-j-r+1),min(param.P,int(p+floor(rho_))),t+1));
                    }
                }

                fTerm[0] = fTerm[0] + act[u-1]*mb_*(1-h(s,p,t))*sum1;
                fTerm[1] = fTerm[1] + act[u-1]*mb_*h(s,p,t)*sum2;
                fTerm[2] = fTerm[2] + act[u-1]*mb_*(1-h(s,p,t))*sum3;
                fTerm[3] = fTerm[3] + act[u-1]*mb_*h(s,p,t)*sum4;
                fTerm[4] = fTerm[4] + act[u-1]*(1-mb_)*(1-h(s,p,t))*sum5;
                fTerm[5] = fTerm[5] + act[u-1]*(1-mb_)*h(s,p,t)*sum6;
            }

            fTerm[6] = g2(s,p,t)*(1-param.mh2)*vf[t+1];
            //for (int i=0; i<7;i++) cout << "fterm["<< i<<"]: " << fTerm[i] << " ";
            double SUM=0;
            for (int i = 0; i<7; i++)
            {
                SUM = SUM + fTerm[i];
            }
            vh(s,p,t) = SUM;
            //cout << endl << "RV of helpers t= " << t << " s= " << s << " p= " << p << "vh="<< vh(s,p,t) << endl;
        }
    }


}

void calcHelpersFoc(const Data& param, Array3& h,Array3& g2,int choice,double* vf, Array3& vhi, Array3& vh, Array3& vb,Array3& mb,Array3& mh,Array3& rho, int t)
{

    double mh_foc=0;
    double g2_foc=0;

    double fTerm[7]; //terms in the helper reproductive value equation
    memset(fTerm, 0, 7*sizeof(double)); //initialize all 0
    double sum1, sum2, sum3, sum4, sum5, sum6;

    Array2 phi(param.S+1,param.P+1); //nest progression probability
    for(int s = 0; s < param.S + 1; s++)
    {
        for(int p = 0; p < param.P + 1; p++)
            phi(s,p)=rho(s,p,t)-floor(rho(s,p,t));
    }



    for (int s = 2; s < param.S + 1; s++) //there are no helpers in colonies of size 0 and 1 (so they are set to 0 by default)
    {
        for(int p = 0; p < param.P + 1; p++)
        {

            sum1=0;
            sum2=0;
            sum3=0;
            sum4=0;
            sum5=0;
            sum6=0;
            for(int j = 0; j < s-1; j++)
            {
                for(int r = 0; r < s-j-1; r++)
                {
                    double ph = bin(s-2,j,mh(s,p,t))*bin(s-2-j,r,g2(s,p,t));
                    sum1 = sum1 + ph * (1.0/double(s-j-r-1))             * (phi(s,p)*vb(s-j-r-1,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)              + (1-phi(s,p))*vb(s-j-r-1,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum2 = sum2 + ph * (1.0/double(s-j-r))               * (phi(s,p)*vb(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)                + (1-phi(s,p))*vb(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum3 = sum3 + ph * (double(s-j-r-2)/double(s-j-r-1)) * (phi(s,p)*vh(s-j-r-1,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)              + (1-phi(s,p))*vh(s-j-r-1,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum4 = sum4 + ph * (double(s-j-r-1)/double(s-j-r))   * (phi(s,p)*vh(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)                + (1-phi(s,p))*vh(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum5 = sum5 + ph *                                     (phi(s,p)*vh(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)                + (1-phi(s,p))*vh(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum6 = sum6 + ph *                                     (phi(s,p)*vh(min(param.S,s-j-r+1),min(param.P,int(p+ceil(rho(s,p,t)))),t+1) + (1-phi(s,p))*vh(min(param.S,s-j-r+1),min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                }
            }

            fTerm[0] = mb(s,p,t)*(1-h(s,p,t))*sum1;
            fTerm[1] = mb(s,p,t)*h(s,p,t)*sum2;
            fTerm[2] = mb(s,p,t)*(1-h(s,p,t))*sum3;
            fTerm[3] = mb(s,p,t)*h(s,p,t)*sum4;
            fTerm[4] = (1-mb(s,p,t))*(1-h(s,p,t))*sum5;
            fTerm[5] = (1-mb(s,p,t))*h(s,p,t)*sum6;
            fTerm[6] = vf[t+1];

//calculate RV of a specific subordinate (for decision making calculations)

            if (choice==1) mh_foc=param.mh1;
            else if (choice==2) mh_foc=param.mh2;
            else if (choice==3) mh_foc=param.mh3;
            else if (choice==4) mh_foc=param.mh4;

            g2_foc=(choice==2 ? 1:0);


            vhi(s,p,t) = (1-mh_foc)*((1.0-g2_foc)*(fTerm[0]+fTerm[1]+fTerm[2]+fTerm[3]+fTerm[4]+fTerm[5])+g2_foc*(1-param.mh2)*fTerm[6]);
        }
    }


}

void calcColonies(const Data& param, Array3& h, Array2& zi, Array3& g1,Array3& g2,Array3& g4, Array3& vc,Array3& mb,Array3& mh,Array3& rho,int t)
{

    double fTerm[6]; //terms in the colony reproductive value equation
    memset(fTerm, 0, 6*sizeof(double)); //initialize all 0

    Array2 phi(param.S+1,param.P+1); //nest progression probability
    for(int s = 0; s < param.S + 1; s++)
    {
        for(int p = 0; p < param.P + 1; p++)
            phi(s,p)=rho(s,p,t)-floor(rho(s,p,t));
    }

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

    //nests of size 0 and 1
    for(int p = 0; p < param.P + 1; p++)
    {
        fTerm[0] = h(0,p,t)*vc(1,p,t+1); //floater joins
        fTerm[1] = (1.0-h(0,p,t))*vc(0,p,t+1); //floater does not join
        vc(0,p,t) = fTerm[0]+fTerm[1];

        fTerm[0] = (1.0-zi(p,t))*param.ms1*h(1,p,t)*(phi(1,p)*(vc(1,min(param.P,int(p+ceil(rho(1,p,t)))),t+1)) + (1-phi(1,p))*vc(1,min(param.P,int(p+floor(rho(1,p,t)))),t+1)); //breeder stays, but dies, floater joins
        fTerm[1] = (1.0-zi(p,t))*param.ms1*(1-h(1,p,t))*(phi(1,p)*vc(0,min(param.P,int(p+ceil(rho(1,p,t)))),t+1) + (1-phi(1,p))*vc(0,min(param.P,int(p+floor(rho(1,p,t)))),t+1)); //breeder stays, but dies, floater does not join
        fTerm[2] = zi(p,t)*h(1,p,t)*vc(1,p,t+1); //breeder leaves, floater joins
        fTerm[3] = zi(p,t)*(1-h(1,p,t))*vc(0,p,t+1); //breeder leaves, floater does not join
        fTerm[4] = (1.0-zi(p,t))*(1-param.ms1)*h(1,p,t)*(phi(1,p)*(vc(2,min(param.P,int(p+ceil(rho(1,p,t)))),t+1)) + (1-phi(1,p))*vc(2,min(param.P,int(p+floor(rho(1,p,t)))),t+1));
        fTerm[5] = (1.0-zi(p,t))*(1-param.ms1)*(1-h(1,p,t))*(phi(1,p)*(vc(1,min(param.P,int(p+ceil(rho(1,p,t)))),t+1)) + (1-phi(1,p))*vc(1,min(param.P,int(p+floor(rho(1,p,t)))),t+1));

        vc(1,p,t) = (1.0-zi(p,t))*(p==param.P?f[1]:0) + fTerm[0] + fTerm[1] + fTerm[2] + fTerm[3] + fTerm[4] + fTerm[5];

    }

    memset(fTerm, 0, 6*sizeof(double)); //initialize all 0
    double sum1, sum2, sum3, sum4;

//nests of size 2 and larger
    for (int s = 2; s < param.S + 1; s++)
    {
        for(int p = 0; p < param.P + 1; p++)
        {

            sum1=0;
            sum2=0;
            sum3=0;
            sum4=0;

            for(int j = 0; j < s; j++)
            {
                for(int r = 0; r < s-j; r++)
                {
                    double ph = bin(s-1,j,mh(s,p,t))*bin(s-1-j,r,g2(s,p,t));
                    sum1 = sum1 + ph * (phi(s,p)*vc(min(param.S,s-j-r+1),min(param.P,int(p+ceil(rho(s,p,t)))),t+1) + (1-phi(s,p))*vc(min(param.S,s-j-r+1),min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum2 = sum2 + ph * (phi(s,p)*vc(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)                + (1-phi(s,p))*vc(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum3 = sum3 + ph * (phi(s,p)*vc(s-j-r-1,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)              + (1-phi(s,p))*vc(s-j-r-1,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum4 = sum4 + ph * (phi(s,p)*vc(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)                + (1-phi(s,p))*vc(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                }
            }

            fTerm[0] = (1.0-mb(s,p,t))*h(s,p,t)*sum1;
            fTerm[1] = (1.0-mb(s,p,t))*(1-h(s,p,t))*sum2;
            fTerm[2] = mb(s,p,t)*(1-h(s,p,t))*sum3;
            fTerm[3] = mb(s,p,t)*h(s,p,t)*sum4;

            vc(s,p,t) = (p==param.P?f[s]:0) + fTerm[0] + fTerm[1] + fTerm[2] + fTerm[3];
            //if (t==17 && p==0 && s==2) {
            //cout << "v(s=" << s << ",p=" << p << ",t=" << t << "): " << vc(s,p,t) << endl;
            //cout << (p==param.P?f[s]:0) << " " << fTerm[0]  << " " <<  fTerm[1]  << " " <<  fTerm[2]  << " " <<  fTerm[3]  << " " <<  fTerm[4] << endl;
            //}
        }
    }


}


void decisionsFloaters(const Data& param, double* vf, double* Vf, Array3& pi, Array3& newpi, Array2& zi, Array3& g2, Array3& vb, Array3& Vb, Array3& vh, Array3& Vh,Array3& mb,Array3& mh,Array3& rho, int sim)
{

//Containers for decision making process
    double Vl[param.T]; //inclusive RV of leaving
    memset(Vl, 0, param.T*sizeof(double)); //initialize all 0
    double vl[param.T]; //inclusive RV of leaving
    memset(vl, 0, param.T*sizeof(double)); //initialize all 0

    Array3 Vs(param.S+1,param.P+1,param.T); //inclusive RV of staying at x of size s and phase p at time t
    Array3 vs(param.S+1,param.P+1,param.T); //inclusive RV of staying at x of size s and phase p at time t

    Array3 VBest(param.S+1,param.P+1,param.T); //inclusive RV of the best decision for floaters

    Array3 phi(param.S+1,param.P+1,param.T); //nest progression probability

    for (int t = 0; t < param.T; t++)
    {
        for(int s = 0; s < param.S + 1; s++)
        {
            for(int p = 0; p < param.P + 1; p++)
                phi(s,p,t)=rho(s,p,t)-floor(rho(s,p,t));
        }
    }

    double term[9]; //terms in the floater reproductive value equation
    memset(term, 0, 9*sizeof(double)); //initialize all 0

//RV of last time step
    Vl[param.T-1] = 0;
    vl[param.T-1] = 0;
    for(int s = 0; s < param.S; s++)
    {
        for(int p = 0; p < param.P + 1; p++)
        {
            Vs(s,p,param.T-1)=0;
            vs(s,p,param.T-1)=0;
        }

    }

    //Nests of size 0 and 1
    for (int t = 0; t < param.T-1; t++)
    {
        Vl[t] = Vf[t+1]; //Reproductive value of leaving is equal with the RV of y for the next time step
        vl[t] = vf[t+1]; //Reproductive value of leaving is equal with the RV of y for the next time step

        for (int p = 0; p < param.P + 1; p++)
        {

            Vs(0,p,t) = Vb(1,p,t+1); //Staying at the empty x
            vs(0,p,t) = vb(1,p,t+1); //Staying at the empty x


            term[0] = (1-zi(p,t))*param.ms1     * (phi(1,p,t)*vb(1,min(param.P,int(p+ceil(rho(1,p,t)))),t+1)   + (1-phi(1,p,t))*vb(1,min(param.P,int(p+floor(rho(1,p,t)))),t+1)); //breeders decides to stay but dies
            term[1] = (1-zi(p,t))*(1-param.ms1) * (phi(1,p,t)*vh(2,min(param.P,int(p+ceil(rho(1,p,t)))),t+1)   + (1-phi(1,p,t))*vh(2,min(param.P,int(p+floor(rho(1,p,t)))),t+1)); //breeder decides to stay and survives
            term[2] = zi(p,t)*vb(1,p,t+1); //breeder decides to leave
            vs(1,p,t) = term[0] + term[1] + term[2];

            term[0] = (1-zi(p,t))*param.ms1     * (phi(1,p,t)*Vb(1,min(param.P,int(p+ceil(rho(1,p,t)))),t+1)   + (1-phi(1,p,t))*Vb(1,min(param.P,int(p+floor(rho(1,p,t)))),t+1)); //breeders decides to stay but dies
            term[1] = (1-zi(p,t))*(1-param.ms1) * (phi(1,p,t)*Vh(2,min(param.P,int(p+ceil(rho(1,p,t)))),t+1)   + (1-phi(1,p,t))*Vh(2,min(param.P,int(p+floor(rho(1,p,t)))),t+1)); //breeder decides to stay and survives
            term[2] = zi(p,t)*Vb(1,p,t+1); //breeder decides to leave
            Vs(1,p,t) = term[0] + term[1] + term[2];

        }
    }


    //Nest of size 2 and higher
    for (int t = 0; t < param.T-1; t++)
    {
        for (int s = 2; s < param.S; s++)    //goes up to S-1, since y can not stay in xs of maximum size S
        {
            for (int p = 0; p < param.P + 1; p++)
            {


                double sum1 = 0;
                double sum2 = 0;
                for(int j = 0; j < s; j++)
                {
                    for(int r = 0; r < s-j; r++)
                    {
                        double ph = bin(s-1,j,mh(s,p,t))*bin(s-1-j,r,g2(s,p,t));
                        term[3] = ph * (1.0/double(s-j-r))             * (phi(s,p,t)*Vb(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)   + (1-phi(s,p,t))*Vb(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1)); //Staying and becoming a breeder
                        term[4] = ph * (double(s-j-r-1)/double(s-j-r)) * (phi(s,p,t)*Vh(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)   + (1-phi(s,p,t))*Vh(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1)); //Staying and becoming a helper, current breeder dies
                        term[5] = ph *                                   (phi(s,p,t)*Vh(s-j-r+1,min(param.P,int(p+ceil(rho(s,p,t)))),t+1) + (1-phi(s,p,t))*Vh(s-j-r+1,min(param.P,int(p+floor(rho(s,p,t)))),t+1)); //Staying and becoming a helper, current breeder does not die
                        term[6] = ph * (1.0/double(s-j-r))             * (phi(s,p,t)*vb(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)   + (1-phi(s,p,t))*vb(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1)); //Staying and becoming a breeder
                        term[7] = ph * (double(s-j-r-1)/double(s-j-r)) * (phi(s,p,t)*vh(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)   + (1-phi(s,p,t))*vh(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1)); //Staying and becoming a helper, current breeder dies
                        term[8] = ph *                                   (phi(s,p,t)*vh(s-j-r+1,min(param.P,int(p+ceil(rho(s,p,t)))),t+1) + (1-phi(s,p,t))*vh(s-j-r+1,min(param.P,int(p+floor(rho(s,p,t)))),t+1)); //Staying and becoming a helper, current breeder does not die
                        //cout << "t: " << t<< " s: " << s << " p: " << p << " j: " << j << " t2: " << term2 << " t3: " << term3 << " t4: " << term4 << endl;
                        sum1 = sum1 + mb(s,p,t)*(term[3]+term[4]) + (1.0-mb(s,p,t))*term[5];
                        sum2 = sum2 + mb(s,p,t)*(term[6]+term[7]) + (1.0-mb(s,p,t))*term[8];
                    }
                }
                Vs(s,p,t) = sum1;
                vs(s,p,t) = sum2;

            }
        }
    }

    //Averaging over p
    Array2 vs_s(param.S+1,param.T);
    Array2 Vs_s(param.S+1,param.T);





    if (sim == param.N-1)
    {

        // Averaging over p
        double zz1, zz2;
        for(int t=0; t<param.T; t++)
        {
            for(int s=0; s<param.S+1; s++)
            {
                zz1=0;
                zz2=0;
                for(int p=0; p<param.P+1; p++)
                {
                    zz1 = zz1 + vs(s,p,t);
                    zz2 = zz2 + Vs(s,p,t);
                }
                vs_s(s,t)=zz1/double(param.P+1);
                Vs_s(s,t)=zz2/double(param.P+1);
            }
        }

        ofstream fitness;
        fitness.open ("fitness.txt", ios::app);
        for(int t=0; t<param.T; t++)
        {
            for(int p=0; p<param.P+1; p++)
            {
                for(int s=0; s<param.S+1; s++)
                    fitness << vs(s,p,t) << " ";
            }
        }
        fitness << endl;
        for(int t=0; t<param.T; t++)
        {
            for(int p=0; p<param.P+1; p++)
            {
                for(int s=0; s<param.S+1; s++)
                    fitness << Vs(s,p,t)-vs(s,p,t) << " ";
            }
        }
        fitness << endl;
        fitness.close();

        ofstream fitnessAverage;
        fitnessAverage.open ("fitnessAverage.txt", ios::app);
        for(int t=0; t<param.T-1; t++)
        {
            for(int s=0; s<param.S+1; s++)
                fitnessAverage << vs_s(s,t) << " ";
        }
        fitnessAverage << endl;
        for(int t=0; t<param.T-1; t++)
        {
            for(int s=0; s<param.S+1; s++)
                fitnessAverage << Vs_s(s,t)-vs_s(s,t) << " ";
        }
        fitnessAverage << endl;
        fitnessAverage.close();

    ofstream ratioFitness;
    double fitRatio=0;
    ratioFitness.open ("ratio.txt", ios::app);
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++) {
    fitRatio = (pi(s,p,t)*vs(s,p,t)+(1-pi(s,p,t))*vl[t])/(pi(s,p,t)*Vs(s,p,t)+(1-pi(s,p,t))*Vl[t]);
    ratioFitness  << fitRatio << " ";
    }}}
    ratioFitness << endl;
    }


//Find the best pi
    for (int t = 0; t < param.T; t++)
    {
        for (int s = 0; s < param.S; s++)
        {
            for (int p = 0; p < param.P + 1; p++)
            {
                VBest(s,p,t) = max(Vs(s,p,t),Vl[t]);
            }
        }
    }



    for (int t = 0; t < param.T; t++)
    {
        for (int s = 0; s < param.S; s++)
        {
            for (int p = 0; p < param.P + 1; p++)
            {
                double L1 = exp(-1.0*(VBest(s,p,t)-Vs(s,p,t))/param.delta);
                double L2 = exp(-1.0*(VBest(s,p,t)-Vl[t])/param.delta);
                newpi(s,p,t) = L1/(L1+L2);
                newpi(param.S,p,t) = 0; //females will not join xs with maximum size
            }
        }
    }

    for (int t = 0; t < param.T; t++)
    {
        for (int s = 0; s < param.S; s++)
        {
            for (int p = 0; p < param.P + 1; p++)
            {
                for (int av=0; av<(sim>param.cN?param.convAveraging:1); av++ ) {

                newpi(s,p,t) = (pi(s,p,t)+newpi(s,p,t))/2.0;

                }


            }

        }
    }

//cout << "pi" << endl;
//printArray3(param,newpi,param.T-1);

}

void decisionsBreeders(const Data& param, double* Vf, Array3& pi, Array3& newpi, Array2& zi, Array3& g2, Array3& vb, Array3& Vb, Array3& vh, Array3& Vh,Array3& mb,Array3& mh,Array3& rho, int sim)
{




}


void decisionsHelpers(const Data& param, Array3& h, double* vf, Array3& nest, Array3& pi,  Array3& g1, Array3& g2, Array3& g3, Array3& g4, Array3& newg1, Array3& newg2, Array3& newg3, Array3& newg4, Array3& vb, Array3& vh, Array3& vc,Array3& mb,Array3& mh,Array3& rho, int sim)
{

//inclusive RV of subordinates choosing the action u1, u2, u3, u4
    Array3 Vh1(param.S+1,param.P+1,param.T);
    Array3 Vh2(param.S+1,param.P+1,param.T);
    Array3 Vh3(param.S+1,param.P+1,param.T);
    Array3 Vh4(param.S+1,param.P+1,param.T);

    Array3 VBest(param.S+1,param.P+1,param.T);


//direct RV of subordinates choosing the action u1, u2, u3, u4
    Array3 vh1(param.S+1,param.P+1,param.T);
    Array3 vh2(param.S+1,param.P+1,param.T);
    Array3 vh3(param.S+1,param.P+1,param.T);
    Array3 vh4(param.S+1,param.P+1,param.T);

//direct RV of colonies with subordinates choosing the action u1, u2, u3, u4
    Array3 vc1(param.S+1,param.P+1,param.T);
    Array3 vc2(param.S+1,param.P+1,param.T);
    Array3 vc3(param.S+1,param.P+1,param.T);
    Array3 vc4(param.S+1,param.P+1,param.T);

//adjusted values of g1
    Array3 g1_foc1(param.S+1,param.P+1,param.T);
    Array3 g1_foc2(param.S+1,param.P+1,param.T);
    Array3 g1_foc3(param.S+1,param.P+1,param.T);
    Array3 g1_foc4(param.S+1,param.P+1,param.T);

//adjusted values of mb
    Array3 mb_foc1(param.S+1,param.P+1,param.T);
    Array3 mb_foc2(param.S+1,param.P+1,param.T);
    Array3 mb_foc3(param.S+1,param.P+1,param.T);
    Array3 mb_foc4(param.S+1,param.P+1,param.T);


    Array3 rho_foc1(param.S+1,param.P+1,param.T);
    Array3 rho_foc2(param.S+1,param.P+1,param.T);
    Array3 rho_foc3(param.S+1,param.P+1,param.T);
    Array3 rho_foc4(param.S+1,param.P+1,param.T);
//There are 4 choices for the focal individual
//u1-usurp,u2-leave,u3-sit and wait,u4-forage



    for (int t = 0; t < param.T; t++)
    {
        for(int s = 2; s < param.S+1; s++)
        {
            for(int p = 0; p < param.P+1; p++)
            {
                g1_foc1(s,p,t) = g1(s,p,t)*double(s-2)/double(s-1) + 1.0/double(s-1);
                g1_foc2(s,p,t) = g1(s,p,t)*double(s-2)/double(s-1);
                g1_foc3(s,p,t) = g1(s,p,t)*double(s-2)/double(s-1);
                g1_foc4(s,p,t) = g1(s,p,t)*double(s-2)/double(s-1);

                if (param.agg==0)
                {
                    mb_foc1(s,p,t) = param.mbi + param.km*g1_foc1(s,p,t);
                    mb_foc2(s,p,t) = param.mbi + param.km*g1_foc2(s,p,t);
                    mb_foc3(s,p,t) = param.mbi + param.km*g1_foc3(s,p,t);
                    mb_foc4(s,p,t) = param.mbi + param.km*g1_foc4(s,p,t);

                    rho_foc1(s,p,t) = param.p0*(1 - param.k*g1_foc1(s,p,t) + param.bp*(g4(s,p,t)*(s-2)));
                    rho_foc2(s,p,t) = param.p0*(1 - param.k*g1_foc2(s,p,t) + param.bp*(g4(s,p,t)*(s-2)));
                    rho_foc3(s,p,t) = param.p0*(1 - param.k*g1_foc3(s,p,t) + param.bp*(g4(s,p,t)*(s-2)));
                    rho_foc4(s,p,t) = param.p0*(1 - param.k*g1_foc4(s,p,t) + param.bp*(g4(s,p,t)*(s-2)+1));
                }

                else
                {

                    mb_foc1(s,p,t) = param.mbi + param.km*(g1(s,p,t)*double(s-2)+1.0)/double(param.S-1);
                    mb_foc2(s,p,t) = param.mbi + param.km*g1(s,p,t)*double(s-2)/double(param.S-1);
                    mb_foc3(s,p,t) = param.mbi + param.km*g1(s,p,t)*double(s-2)/double(param.S-1);
                    mb_foc4(s,p,t) = param.mbi + param.km*g1(s,p,t)*double(s-2)/double(param.S-1);

                    rho_foc1(s,p,t) = param.p0*(1-param.k*(g1(s,p,t)*double(s-2)+1.0)/double(param.S-1) + param.bp*(g4(s,p,t)*(s-2)));
                    rho_foc2(s,p,t) = param.p0*(1-param.k*g1(s,p,t)*double(s-2)/double(param.S-1)       + param.bp*(g4(s,p,t)*(s-2)));
                    rho_foc3(s,p,t) = param.p0*(1-param.k*g1(s,p,t)*double(s-2)/double(param.S-1)       + param.bp*(g4(s,p,t)*(s-2)));
                    rho_foc4(s,p,t) = param.p0*(1-param.k*g1(s,p,t)*double(s-2)/double(param.S-1)       + param.bp*(g4(s,p,t)*(s-2)+1));



                }


            }
        }
    }


//Last time T-1
    for(int s = 2; s < param.S+1; s++)
    {
        for(int p = 0; p < param.P+1; p++)
        {
            vh1(s,p,param.T-1) = 0;
            vh2(s,p,param.T-1) = 0;
            vh3(s,p,param.T-1) = 0;
            vh4(s,p,param.T-1) = 0;
        }
    }

    for(int s = 1; s < param.S+1; s++)
    {
        for(int p = 0; p < param.P+1; p++)
        {


            if (s>1)
            {
                if (param.agg==0)
                {
                    vc1(s,p,param.T-1) = (p==param.P?param.f0*(1.0-param.k*g1_foc1(s,p,param.T-1)+param.bf*(g4(s,p,param.T-1)*(s-2))):0);
                    vc2(s,p,param.T-1) = (p==param.P?param.f0*(1.0-param.k*g1_foc2(s,p,param.T-1)+param.bf*(g4(s,p,param.T-1)*(s-2))):0);
                    vc3(s,p,param.T-1) = (p==param.P?param.f0*(1.0-param.k*g1_foc3(s,p,param.T-1)+param.bf*(g4(s,p,param.T-1)*(s-2))):0);
                    vc4(s,p,param.T-1) = (p==param.P?param.f0*(1.0-param.k*g1_foc4(s,p,param.T-1)+param.bf*(g4(s,p,param.T-1)*(s-2)+1)):0);
                }
                else
                {
                    vc1(s,p,param.T-1) = (p==param.P?param.f0*(1.0-param.k*(g1(s,p,param.T-1)*double(s-2)+1.0)/double(param.S-1) + param.bf*(g4(s,p,param.T-1)*(s-2))):0);
                    vc2(s,p,param.T-1) = (p==param.P?param.f0*(1.0-param.k*g1(s,p,param.T-1)*double(s-2)/double(param.S-1)       + param.bf*(g4(s,p,param.T-1)*(s-2))):0);
                    vc3(s,p,param.T-1) = (p==param.P?param.f0*(1.0-param.k*g1(s,p,param.T-1)*double(s-2)/double(param.S-1)       + param.bf*(g4(s,p,param.T-1)*(s-2))):0);
                    vc4(s,p,param.T-1) = (p==param.P?param.f0*(1.0-param.k*g1(s,p,param.T-1)*double(s-2)/double(param.S-1)       + param.bf*(g4(s,p,param.T-1)*(s-2)+1)):0);
                }
            }
            else
            {
                vc1(s,p,param.T-1) = (p==param.P?param.f0:0);
                vc2(s,p,param.T-1) = (p==param.P?param.f0:0);
                vc3(s,p,param.T-1) = (p==param.P?param.f0:0);
                vc4(s,p,param.T-1) = (p==param.P?param.f0:0);
            }

        }
    }




// Time steps 0 .. T-2
    for (int t=param.T-2; t>=0; t--)
    {

        calcHelpersFoc(param,h,g2,1,vf,vh1,vh,vb,mb_foc1,mh,rho_foc1,t); //calculates direct RV of action 1 vh1
        calcHelpersFoc(param,h,g2,2,vf,vh2,vh,vb,mb_foc2,mh,rho_foc2,t); //calculates direct RV of action 2 vh2
        calcHelpersFoc(param,h,g2,3,vf,vh3,vh,vb,mb_foc3,mh,rho_foc3,t); //calculates direct RV of action 3 vh3
        calcHelpersFoc(param,h,g2,4,vf,vh4,vh,vb,mb_foc4,mh,rho_foc4,t); //calculates direct RV of action 4 vh4

        calcColDecision(param,h,g1,1,g2,g4,vc1,vc,mb_foc1,mh,rho_foc1,t); //calculates the direct RV of a colony that contains the focal individual choosing the action 1 vc1
        calcColDecision(param,h,g1,2,g2,g4,vc2,vc,mb_foc2,mh,rho_foc2,t); //calculates the direct RV of a colony that contains the focal individual choosing the action 2 vc1
        calcColDecision(param,h,g1,3,g2,g4,vc3,vc,mb_foc3,mh,rho_foc3,t); //calculates the direct RV of a colony that contains the focal individual choosing the action 3 vc1
        calcColDecision(param,h,g1,4,g2,g4,vc4,vc,mb_foc4,mh,rho_foc4,t); //calculates the direct RV of a colony that contains the focal individual choosing the action 4 vc1

    }



    for (int t = 0; t < param.T; t++)
    {
        for(int s = 2; s < param.S + 1; s++)
        {
            for(int p = 0; p < param.P + 1; p++)
            {
                Vh1(s,p,t) = vh1(s,p,t) + (vc1(s,p,t)-vc(s-1,p,t)-vh1(s,p,t))*param.r;
                Vh2(s,p,t) = vh2(s,p,t) + (vc2(s,p,t)-vc(s-1,p,t)-vh2(s,p,t))*param.r;
                Vh3(s,p,t) = vh3(s,p,t) + (vc3(s,p,t)-vc(s-1,p,t)-vh3(s,p,t))*param.r;
                Vh4(s,p,t) = vh4(s,p,t) + (vc4(s,p,t)-vc(s-1,p,t)-vh4(s,p,t))*param.r;
            }
        }
    }




    /*
        for (int t=0;t<18;t++) {

        cout << "vc(t="<< t<<")" << endl;
        printArray3(param,vc,t);

        cout << "vc1(t="<< t<<")" << endl;
        printArray3(param,vc1,t);
        cout << "vc2(t="<< t<<")" << endl;
        printArray3(param,vc2,t);
        cout << "vc3(t="<< t<<")" << endl;
        printArray3(param,vc3,t);
        cout << "vc4(t="<< t<<")" << endl;
        printArray3(param,vc4,t);


        cout << "vh1(t="<< t<<")" << endl;
        printArray3(param,vh1,t);
        cout << "vh2(t="<< t<<")" << endl;
        printArray3(param,vh2,t);
        cout << "vh3(t="<< t<<")" << endl;
        printArray3(param,vh3,t);
        cout << "vh4(t="<< t<<")" << endl;
        printArray3(param,vh4,t);

        cout << "Vh1(t="<< t<<")" << endl;
        printArray3(param,Vh1,t);
        cout << "Vh2(t="<< t<<")" << endl;
        printArray3(param,Vh2,t);
        cout << "Vh3(t="<< t<<")" << endl;
        printArray3(param,Vh3,t);
        cout << "Vh4(t="<< t<<")" << endl;
        printArray3(param,Vh4,t);



        }
    */








    double max1, max2;

//Find the best strategy
    for (int t = 0; t < param.T; t++)
    {
        for (int s = 2; s < param.S + 1; s++)
        {
            for (int p = 0; p < param.P + 1; p++)
            {
                max1 = max(Vh1(s,p,t),Vh2(s,p,t));
                max2 = max(Vh3(s,p,t),Vh4(s,p,t));
                VBest(s,p,t) = max(max1,max2);
            }
        }
    }

    double gSum = 0;

    for (int t = 0; t < param.T; t++)
    {
        for (int s = 2; s < param.S + 1; s++)
        {
            for (int p = 0; p < param.P + 1; p++)
            {
                double L1 = exp(-1.0*(VBest(s,p,t)-Vh1(s,p,t))/param.delta);
                double L2 = exp(-1.0*(VBest(s,p,t)-Vh2(s,p,t))/param.delta);
                double L3 = exp(-1.0*(VBest(s,p,t)-Vh3(s,p,t))/param.delta);
                double L4 = exp(-1.0*(VBest(s,p,t)-Vh4(s,p,t))/param.delta);
                newg1(s,p,t) = L1/(L1+L2+L3+L4);
                newg2(s,p,t) = L2/(L1+L2+L3+L4);
                newg3(s,p,t) = L3/(L1+L2+L3+L4);
                newg4(s,p,t) = L4/(L1+L2+L3+L4);

                /*
                if(abs(g1(s,p,t)-newg1(s,p,t)) > param.convPrec || abs(g2(s,p,t)-newg2(s,p,t)) > param.convPrec || abs(g3(s,p,t)-newg3(s,p,t)) > param.convPrec ||abs(g4(s,p,t)-newg4(s,p,t)) > param.convPrec ) {

                cout << "t=" << t << " s=" << s << " p=" << p << endl;
                cout << "dif " << abs(g1(s,p,t)-newg1(s,p,t)) << " "
                     << abs(g2(s,p,t)-newg2(s,p,t)) << " "
                     << abs(g3(s,p,t)-newg3(s,p,t)) << " "
                     << abs(g4(s,p,t)-newg4(s,p,t)) << endl;
                cout << " g "<< g1(s,p,t) << " " << g2(s,p,t) << " " << g3(s,p,t) << " " << g4(s,p,t) << " " << endl;
                cout << "newg "<< newg1(s,p,t) << " " << newg2(s,p,t) << " " << newg3(s,p,t) << " " << newg4(s,p,t) << " " << endl;
                cout << "L1=" << L1 << " L2=" << L2 << " L3=" << L3 << " L4=" << L4 << endl;
                cout << "VH1=" << Vh1(s,p,t) << " VH2=" << Vh2(s,p,t) << " VH3=" << Vh3(s,p,t) << " VH4=" << Vh4(s,p,t) << endl;

                }
                */

                gSum = newg1(s,p,t) + newg2(s,p,t) + newg3(s,p,t) + newg4(s,p,t);
                if (!almostEqual(gSum, 1, 0.000000001)) cout << "Error: Subordinate choices do not add up to 1! gSum (t="<< t << ")= " << gSum << endl;
            }
        }
    }

    for (int t = 0; t < param.T; t++)
    {
        for (int s = 0; s < 2; s++)
        {
            for (int p = 0; p < param.P + 1; p++)
            {
                newg1(s,p,t) = 0;
                newg2(s,p,t) = 0;
                newg3(s,p,t) = 0;
                newg4(s,p,t) = 0;
            }
        }
    }

    for (int t = 0; t < param.T; t++)
    {
        for (int s = 2; s < param.S + 1; s++)
        {
            for (int p = 0; p < param.P + 1; p++)
            {


                for (int av=0; av<(sim>param.cN?param.convAveraging:1); av++ ) {

                newg1(s,p,t) = (g1(s,p,t)+newg1(s,p,t))/2.0;
                newg2(s,p,t) = (g2(s,p,t)+newg2(s,p,t))/2.0;
                newg3(s,p,t) = (g3(s,p,t)+newg3(s,p,t))/2.0;
                newg4(s,p,t) = (g4(s,p,t)+newg4(s,p,t))/2.0;

                }



            }

        }
    }



//for (int t = 0; t < param.T; t++) {
//cout << "g1: " << t << endl;
//printArray3(param,newg1,t);
//cout << "g2:" << t << endl;
//printArray3(param,newg2,t);
//cout << "g3:" << t << endl;
//printArray3(param,newg3,t);
//cout << "g4:" << t << endl;
//printArray3(param,newg4,t);
//}

    if (sim == param.N-1)
    {
        printFitness(param,vh1,vh2,vh3,vh4,Vh1,Vh2,Vh3,Vh4,g1,g2,g3,g4);
    }

}

void calcColDecision(const Data& param, Array3& h,  Array3& g1, int choice, Array3& g2, Array3& g4, Array3& vci, Array3& vc, Array3& mb, Array3& mh, Array3& rho,int t)
{

    double g2_foc, mh_foc;

    if (choice==1) mh_foc=param.mh1;
    else if (choice==2) mh_foc=param.mh2;
    else if (choice==3) mh_foc=param.mh3;
    else if (choice==4) mh_foc=param.mh4;

    g2_foc=(choice==2 ? 1:0);

    double fTerm[8]; //terms in the colony reproductive value equation

    Array2 phi(param.S+1,param.P+1); //nest progression probability
    for(int s = 0; s < param.S + 1; s++)
    {
        for(int p = 0; p < param.P + 1; p++)
            phi(s,p)=rho(s,p,t)-floor(rho(s,p,t));
    }

    double f[param.S+1]; //reproduction function

    f[0] = 0; //no reproduction in empty nests
    f[1] = param.f0;
    if (param.S > 1)
    {
        for(int s = 2; s < param.S + 1; s++)
        {
            if (param.agg==0) f[s] = param.f0*((1.0-param.k*(g1(s,param.P,t)*double(s-2)/double(s-1) + (choice==1?1.0/double(s-1):0)))+param.bf*(g4(s,param.P,t)*(s-2)+(choice==4?1:0)));
            else f[s] = param.f0*((1.0-param.k*(g1(s,param.P,t)*double(s-2)+(choice==1?1.0:0.0))/double(param.S-1)) + param.bf*(g4(s,param.P,t)*(s-2)+(choice==4?1:0)));
        }
    }


    memset(fTerm, 0, 8*sizeof(double)); //initialize all 0
    double sum1, sum2, sum3, sum4;

    double ph, q1, q2;
    q1 = (1-mh_foc)*(1-g2_foc);
    q2 = mh_foc+(1-mh_foc)*g2_foc;

    //colonies of size 2 and larger, since I do not need to calculate vci for s<2
    for (int s = 2; s < param.S + 1; s++)
    {
        for(int p = 0; p < param.P + 1; p++)
        {

            sum1=0;
            sum2=0;
            sum3=0;
            sum4=0;

            for(int j = 0; j < s-1; j++)
            {
                for(int r = 0; r < s-j-1; r++)
                {
                    ph = bin(s-2,j,mh(s,p,t))*bin(s-2-j,r,g2(s,p,t));
                    sum1 = sum1 + ph * (phi(s,p)*vc(min(param.S,s-j-r+1),min(param.P,int(p+ceil(rho(s,p,t)))),t+1) + (1-phi(s,p))*vc(min(param.S,s-j-r+1),min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum2 = sum2 + ph * (phi(s,p)*vc(s-j-r,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)                + (1-phi(s,p))*vc(s-j-r,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum3 = sum3 + ph * (phi(s,p)*vc(s-j-r-1,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)              + (1-phi(s,p))*vc(s-j-r-1,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                    sum4 = sum4 + ph * (phi(s,p)*vc(s-j-r-2,min(param.P,int(p+ceil(rho(s,p,t)))),t+1)              + (1-phi(s,p))*vc(s-j-r-2,min(param.P,int(p+floor(rho(s,p,t)))),t+1));
                }
            }


            fTerm[0] = q1*(1.0-mb(s,p,t))*h(s,p,t)*sum1;
            fTerm[1] = q2*(1.0-mb(s,p,t))*h(s,p,t)*sum2;
            fTerm[2] = q1*(1.0-mb(s,p,t))*(1-h(s,p,t))*sum2;
            fTerm[3] = q2*(1.0-mb(s,p,t))*(1-h(s,p,t))*sum3;
            fTerm[4] = q1*mb(s,p,t)*(1-h(s,p,t))*sum3;
            fTerm[5] = q2*mb(s,p,t)*(1-h(s,p,t))*sum4;
            fTerm[6] = q1*mb(s,p,t)*h(s,p,t)*sum2;
            fTerm[7] = q2*mb(s,p,t)*h(s,p,t)*sum3;


            vci(s,p,t) = (p==param.P?f[s]:0) + fTerm[0] + fTerm[1] + fTerm[2] + fTerm[3] + fTerm[4] + fTerm[5] + fTerm[6] + fTerm[7];
            //if (t==17 && p==0 &&s==2) {
            //cout << "vci(s=" << s << ",p=" << p << ",t=" << t << "): " << vci(s,p,t) << endl;
            //cout << (p==param.P?f[s]:0) << " " << fTerm[0]  << " " <<  fTerm[1]  << " " <<  fTerm[2]  << " " <<  fTerm[3]  << " " <<  fTerm[4]  << " " <<  fTerm[5]  << " " <<  fTerm[6]  << " " <<  fTerm[7] << endl;
            //}
        }
    }


}

































