#include <iostream>
#include <fstream>
#include "print.h"

using namespace std;


void printParam(const Data& param) {

    ofstream sfile; //delete the content of previously existing data, and write the important parameter values
    sfile.open ("pi.txt", ios::trunc);
    sfile << param.S+1 << " " << param.P+1 << " " << param.T << endl;
    sfile << param.r << " " << param.delta << " " << param.pF << " " << param.N << " " << param.mbi << " " << param.mh1 << " " << param.mh2 << " " << param.mh3 << " " << param.mh4 << " ";
    sfile << param.mf << " " << param.X << " " << param.y0 << " " << param.f0 << " " << param.p0 << " " << param.bf << " ";
    sfile << param.k << " " << param.km << " " << endl;
    sfile.close();

    ofstream g1file; //delete the content of previously existing data, and write the important parameter values
    g1file.open ("g1.txt", ios::trunc);
    g1file << param.S+1 << " " << param.P+1 << " " << param.T << endl;
    g1file << param.r << " " << param.delta << " " << param.pF << " " << param.N << " " << param.mbi << " " << param.mh1 << " " << param.mh2 << " " << param.mh3 << " " << param.mh4 << " ";
    g1file << param.mf << " " << param.X << " " << param.y0 << " " << param.f0 << " " << param.p0 << " " << param.bf << " ";
    g1file << param.k << " " << param.km << " " << endl;
    g1file.close();

    ofstream g2file; //delete the content of previously existing data, and write the important parameter values
    g2file.open ("g2.txt", ios::trunc);
    g2file << param.S+1 << " " << param.P+1 << " " << param.T << endl;
    g2file << param.r << " " << param.delta << " " << param.pF << " " << param.N << " " << param.mbi << " " << param.mh1 << " " << param.mh2 << " " << param.mh3 << " " << param.mh4 << " ";
    g2file << param.mf << " " << param.X << " " << param.y0 << " " << param.f0 << " " << param.p0 << " " << param.bf << " ";
    g2file << param.k << " " << param.km << " " << endl;
    g2file.close();

    ofstream g3file; //delete the content of previously existing data, and write the important parameter values
    g3file.open ("g3.txt", ios::trunc);
    g3file << param.S+1 << " " << param.P+1 << " " << param.T << endl;
    g3file << param.r << " " << param.delta << " " << param.pF << " " << param.N << " " << param.mbi << " " << param.mh1 << " " << param.mh2 << " " << param.mh3 << " " << param.mh4 << " ";
    g3file << param.mf << " " << param.X << " " << param.y0 << " " << param.f0 << " " << param.p0 << " " << param.bf << " ";
    g3file << param.k << " " << param.km << " " << endl;
    g3file.close();

    ofstream g4file; //delete the content of previously existing data, and write the important parameter values
    g4file.open ("g4.txt", ios::trunc);
    g4file << param.S+1 << " " << param.P+1 << " " << param.T << endl;
    g4file << param.r << " " << param.delta << " " << param.pF << " " << param.N << " " << param.mbi << " " << param.mh1 << " " << param.mh2 << " " << param.mh3 << " " << param.mh4 << " ";
    g4file << param.mf << " " << param.X << " " << param.y0 << " " << param.f0 << " " << param.p0 << " " << param.bf << " ";
    g4file << param.k << " " << param.km << " " << endl;
    g4file.close();

    ofstream dfile;
    dfile.open ("forward.txt", ios::trunc);
    dfile << param.S+1 << " " << param.P+1 << " " << param.T  << endl;
    dfile.close();

    ofstream fitness;
    fitness.open ("fitness.txt", ios::trunc);
    fitness << param.S+1 << " " << param.P+1 << " " << param.T << endl;
    fitness.close();

    ofstream fitnessAverage;
    fitnessAverage.open ("fitnessAverage.txt", ios::trunc);
    fitnessAverage << param.S+1 << " " << param.T-1 << endl;
    fitnessAverage.close();

    ofstream ratioFitness;
    fitness.open ("ratio.txt", ios::trunc);
    fitness << param.S+1 << " " << param.P+1 << " " << param.T << endl;
    fitness.close();

    ofstream ratioFitnessAverage;
    fitnessAverage.open ("ratioAverage.txt", ios::trunc);
    fitnessAverage << param.S+1 << " " << param.T-1 << endl;
    fitnessAverage.close();




}

void printData(const Data& param, Array3& nest, Array3& pi, Array3& g1, Array3& g2, Array3& g3, Array3& g4, double* floaters) {

    ofstream strat0;
    strat0.open ("pi.txt", ios::app);
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    strat0 << pi(s,p,t) << " ";
    }}
    strat0.close();

    ofstream strat1;
    strat1.open ("g1.txt", ios::app);
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    strat1 << g1(s,p,t) << " ";
    }}
    strat1.close();

    ofstream strat2;
    strat2.open ("g2.txt", ios::app);
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    strat2 << g2(s,p,t) << " ";
    }}
    strat2.close();

    ofstream strat3;
    strat3.open ("g3.txt", ios::app);
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    strat3 << g3(s,p,t) << " ";
    }}
    strat3.close();

    ofstream strat4;
    strat4.open ("g4.txt", ios::app);
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    strat4 << g4(s,p,t) << " ";
    }}
    strat4.close();



    ofstream dfile;
    dfile.open ("forward.txt", ios::app);
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++){
    dfile << nest(s,p,t) << " ";
    }}}
    dfile << endl;
    for (int t = 0; t < param.T; t++) dfile << floaters[t] << " ";
    dfile << endl;
    dfile.close();

    //Averaging over p
    Array2 pi_s(param.S+1,param.T);
    Array2 g1_s(param.S+1,param.T);
    Array2 g2_s(param.S+1,param.T);
    Array2 g3_s(param.S+1,param.T);
    Array2 g4_s(param.S+1,param.T);

    double sum0, sum1, sum2, sum3, sum4;

    for(int t=0;t<param.T;t++){
    for(int s=0;s<param.S+1;s++){
    sum0=0;sum1=0;sum2=0;sum3=0;sum4=0;
    for(int p=0;p<param.P+1;p++){
    sum0=sum0+pi(s,p,t);
    sum1=sum1+g1(s,p,t);
    sum2=sum2+g2(s,p,t);
    sum3=sum3+g3(s,p,t);
    sum4=sum4+g4(s,p,t);
    }
    pi_s(s,t)=sum0/double(param.P+1);
    g1_s(s,t)=sum1/double(param.P+1);
    g2_s(s,t)=sum2/double(param.P+1);
    g3_s(s,t)=sum3/double(param.P+1);
    g4_s(s,t)=sum4/double(param.P+1);
    }}


    ofstream mean;
    mean.open ("meanS.txt", ios::trunc);
    mean << param.S+1 << " " << param.T << endl;
    mean.close();

    mean.open ("meanS.txt", ios::app);

    for(int t=0;t<param.T;t++){
    for(int s=0;s<param.S;s++)
    mean << pi_s(s,t) << " ";
    }
    mean << endl;
    for(int t=0;t<param.T;t++){
    for(int s=2;s<param.S+1;s++)
    mean << g1_s(s,t) << " ";
    }
    mean << endl;
    for(int t=0;t<param.T;t++){
    for(int s=2;s<param.S+1;s++)
    mean << g2_s(s,t) << " ";
    }
    mean << endl;
    for(int t=0;t<param.T;t++){
    for(int s=2;s<param.S+1;s++)
    mean << g3_s(s,t) << " ";
    }
    mean << endl;
    for(int t=0;t<param.T;t++){
    for(int s=2;s<param.S+1;s++)
    mean << g4_s(s,t) << " ";
    }
    mean << endl;
    mean.close();


}

void printFitness(const Data& param, Array3& vh1, Array3& vh2, Array3& vh3, Array3& vh4, Array3& Vh1, Array3& Vh2, Array3& Vh3, Array3& Vh4, Array3& g1, Array3& g2, Array3& g3, Array3& g4) {

    ofstream fitness;
    fitness.open ("fitness.txt", ios::app);
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    fitness << vh1(s,p,t) << " ";
    }}
    fitness << endl;
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    fitness << Vh1(s,p,t)-vh1(s,p,t) << " ";
    }}
    fitness << endl;

    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    fitness << vh2(s,p,t) << " ";
    }}
    fitness << endl;
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    fitness << Vh2(s,p,t)-vh2(s,p,t) << " ";
    }}
    fitness << endl;

    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    fitness << vh3(s,p,t) << " ";
    }}
    fitness << endl;
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    fitness << Vh3(s,p,t)-vh3(s,p,t) << " ";
    }}
    fitness << endl;

    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    fitness << vh4(s,p,t) << " ";
    }}
    fitness << endl;
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++)
    fitness << Vh4(s,p,t)-vh4(s,p,t) << " ";
    }}
    fitness << endl;
    fitness.close();

    //Averaging over p
    Array2 vh1_s(param.S+1,param.T);
    Array2 Vh1_s(param.S+1,param.T);
    Array2 vh2_s(param.S+1,param.T);
    Array2 Vh2_s(param.S+1,param.T);
    Array2 vh3_s(param.S+1,param.T);
    Array2 Vh3_s(param.S+1,param.T);
    Array2 vh4_s(param.S+1,param.T);
    Array2 Vh4_s(param.S+1,param.T);

    // Averaging over p
    double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8;
    for(int t=0;t<param.T;t++){
    for(int s=0;s<param.S+1;s++){
    sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;sum6=0;sum7=0;sum8=0;
    for(int p=0;p<param.P+1;p++){

    sum1 = sum1 + vh1(s,p,t);
    sum2 = sum2 + Vh1(s,p,t);

    sum3 = sum3 + vh2(s,p,t);
    sum4 = sum4 + Vh2(s,p,t);

    sum5 = sum5 + vh3(s,p,t);
    sum6 = sum6 + Vh3(s,p,t);

    sum7 = sum7 + vh4(s,p,t);
    sum8 = sum8 + Vh4(s,p,t);
    }

    vh1_s(s,t)=sum1/double(param.P+1);
    Vh1_s(s,t)=sum2/double(param.P+1);

    vh2_s(s,t)=sum1/double(param.P+1);
    Vh2_s(s,t)=sum2/double(param.P+1);

    vh3_s(s,t)=sum1/double(param.P+1);
    Vh3_s(s,t)=sum2/double(param.P+1);

    vh4_s(s,t)=sum1/double(param.P+1);
    Vh4_s(s,t)=sum2/double(param.P+1);

    }}

    ofstream fitnessAverage;
    fitnessAverage.open ("fitnessAverage.txt", ios::app);

    for(int t=0;t<param.T-1;t++){
    for(int s=0;s<param.S+1;s++)
    fitnessAverage << vh1_s(s,t) << " ";
    }
    fitnessAverage << endl;
    for(int t=0;t<param.T-1;t++){
    for(int s=0;s<param.S+1;s++)
    fitnessAverage << Vh1_s(s,t)-vh1_s(s,t) << " ";
    }
    fitnessAverage << endl;

    for(int t=0;t<param.T-1;t++){
    for(int s=0;s<param.S+1;s++)
    fitnessAverage << vh2_s(s,t) << " ";
    }
    fitnessAverage << endl;
    for(int t=0;t<param.T-1;t++){
    for(int s=0;s<param.S+1;s++)
    fitnessAverage << Vh2_s(s,t)-vh2_s(s,t) << " ";
    }
    fitnessAverage << endl;

    for(int t=0;t<param.T-1;t++){
    for(int s=0;s<param.S+1;s++)
    fitnessAverage << vh3_s(s,t) << " ";
    }
    fitnessAverage << endl;
    for(int t=0;t<param.T-1;t++){
    for(int s=0;s<param.S+1;s++)
    fitnessAverage << Vh3_s(s,t)-vh3_s(s,t) << " ";
    }
    fitnessAverage << endl;

    for(int t=0;t<param.T-1;t++){
    for(int s=0;s<param.S+1;s++)
    fitnessAverage << vh4_s(s,t) << " ";
    }
    fitnessAverage << endl;
    for(int t=0;t<param.T-1;t++){
    for(int s=0;s<param.S+1;s++)
    fitnessAverage << Vh4_s(s,t)-vh4_s(s,t) << " ";
    }
    fitnessAverage << endl;

    fitnessAverage.close();


    ofstream ratioFitness;
    double fitRatio=0;
    ratioFitness.open ("ratio.txt", ios::app);
    for(int t=0;t<param.T;t++){
    for(int p=0;p<param.P+1;p++){
    for(int s=0;s<param.S+1;s++) {
    fitRatio = (g1(s,p,t)*vh1(s,p,t)+g2(s,p,t)*vh2(s,p,t)+g3(s,p,t)*vh3(s,p,t)+g4(s,p,t)*vh4(s,p,t))/(g1(s,p,t)*Vh1(s,p,t)+g2(s,p,t)*Vh2(s,p,t)+g3(s,p,t)*Vh3(s,p,t)+g4(s,p,t)*Vh4(s,p,t));
    ratioFitness  << fitRatio << " ";
    }}}
    ratioFitness << endl;


    ofstream ratioFitnessAverage;
    ratioFitnessAverage.open ("ratioAverage.txt", ios::app);
    ratioFitnessAverage << param.S+1 << " " << param.T << endl;
    ratioFitnessAverage.close();
}

