// Algorithm obtained from: Villagómez Correa, Stephany Rocío (2016) https://bibdigital.epn.edu.ec/handle/15000/15440 

#include <iostream>
#include <fstream>
#include <stdlib.h>  //required to use 'rand()'
#include <math.h>
#include <time.h>    //required to use 'srand(time(NULL))'
#include <ctime>

using namespace std;

const double dte=0.01;
const int N=500000000;
const double sigm=0.02;
const double chi=24.725;
double pi=3.1415926535897;

// Function that calculates a random number between 1-a/2 and 1+a/2
double randfactor (double a)
{
 double r;
 r= (double) rand()/(double)RAND_MAX;
 r=1+(r-0.5)*a;
 return r;
}

//Reaction time
double rtime (double A, double c)
{
 double r, t;

 r=(double) rand()/RAND_MAX;
 t= (1/(A*c))*log(1/r);

 return t;
}

//Reaction time
double ttime (double A, double c, double B, double U)
{
  double r, t;

  if (B>U)
  {
   t= rtime(A,c);
  }
  else
  {
   t=1000000;
  }
 
  return t;
}

// Total fluorescence
double Fluo (double fo, double a, double b, double c, double d, double e, double f, double f1, double f2, double f3)
{
 double r;
 r=fo+a*f1/(1-b*f1)+c*f2/(1-d*f2)+e*f3/(1-f*f3);
 return r;
}

// First part of the fluorescence curve
double Fluo1 (double fo, double a, double b, double f)
{
 double r;
 r=fo+a*f/(1-b*f);
 return r;
}

// First and second part of the fluorescence curve
double Fluo2 (double fo, double a1, double b1, double a2, double b2, double f1, double f2)
{
 double r;
 r=fo+a1*f1/(1-b1*f1)+a2*f2/(1-b2*f2);
 return r;
}

double Kinetics1 (double& dtt, double& A1, double& B1, double c1)
{
 double t,tao,x;

 t=dtt;
 while( t < dte)
 {
  if (A1>0)
  {
   tao= rtime(A1,c1);
   A1--;
   B1++;
   t=t+tao;
  }
  else
  {
   t=dte;
  }
 }
 dtt=t-dte;
 x=B1/N;
 return x;
}

double Kinetics2 (double& dtt, double B1, double U1, double& A2, double& B2, double c2)
{
 double t,tao,x;

 t=dtt;
 while(t < dte)
 {
  if (A2>0 && B1>U1)
  {
   tao= rtime(A2,c2);
   A2--;
   B2++;
   t=t+tao;
  }
   else
  {
   t=dte;
  }
 }
 dtt=t-dte;

 x=B2/N; 
 return x;
}

void update (int& cn, int& count, double& e, double& ex, double& v, double& vx, double& d, double& dx)
{
 double q=0.5;
 count++;
 if (ex<e)
   {
    e=ex;
    v=vx;
    cn++;
    dx= dx/(1+q*count);
    cout<<count<<" "<<v<<" "<<e<<endl;
   }
   else
   {
    dx= -dx/(1+q*count);
   }
   d=dx;
}

int main()
{
 int U1,U2,count,countx,pts,p,n;
 double nllm, ff, dtt,mm,x1,x2,x3,A1,B1,A2,B2,A3,B3,b1,b2,b3,c1,c2,c3,k1,k2,k3,t,tt,tao,tao1,tao2,tao3,test,a1,aa1,a2,aa2,a3,aa3,p1,pp1,p2,pp2,p3,pp3,F,F1,F2,F3,FF,dif,sq,er,erx,u1,u2,v1,v2,chi2;
 clock_t start_time = clock();

 srand(time(NULL));   //required for "randomness"
 ifstream exp ("data_for_analysis.csv");
 ofstream sim;
 sim.open ("model.dat");
 ofstream par ("parmle.dat");
 ofstream time1 ("time1.dat");
 
 pts=24000; 

 ff=1.7; 

 k1=2.74932;    //First kinetic constant (1/ms)
 a1=0.98108;
 p1=0.46775;

 k2=0.14549;   //Second kinetic constant (1/ms)
 a2=0.92052;
 p2=0.52087;

 k3=0.00936;  //Third kinetic constant (1/ms)
 a3=0.69296;
 p3=0.69845;

 u1=0.94634;
 u2=0.95017;

 double eu,eux,dd,delta;

 par<<"k1     "<<"a1     "<<"p1     "<<"k2     "<<"a2     "<<"p2     "<<"k3     "<<"a3     "<<"p3     "<<"u1   "<<"u2   "<<"nll  "<<endl;

  double dtt1,dtt2,dtt3;
  int cn=0;
  er=0.025;
  dd=0;
  count=0;
  do
  {
   cout<<"Iteration: "<<count+1<<"  ";
   dtt1=0;
   dtt2=0;
   dtt3=0;
   A1=N;
   A2=N;
   A3=N;
   B1=0;
   B2=0;
   B3=0;
   sq=0;

   c1=randfactor(0.01)*k1;
   aa1=randfactor(0.01)*a1;
   pp1=randfactor(0.01)*p1;

   c2=randfactor(0.01)*k2;  
   aa2=randfactor(0.01)*a2;
   pp2=randfactor(0.01)*p2;

   c3=randfactor(0.01)*k3;  
   aa3=randfactor(0.01)*a3;
   pp3=randfactor(0.01)*p3;

   v1=randfactor(0.01)*u1;
   v2=randfactor(0.01)*u2;
   
   U1=v1*N;
   U2=v2*N;

   for( int i = 0; i < pts; i++  )
   {
    b1=Kinetics1 (dtt1, A1, B1, c1);
    b2=Kinetics2 (dtt2, B1, U1, A2, B2, c2);
    b3=Kinetics2 (dtt3, B2, U2, A3, B3, c3);
    exp>>t>>FF;
    F=Fluo (ff, aa1, pp1, aa2, pp2, aa3, pp3, b1, b2, b3);
    dif=FF-F;
    sq=sq+dif*dif;
   }
   chi2=sq/(sigm*sigm);
   erx=sqrt(sq/pts);
   par<<c1<<" "<<aa1<<" "<<pp1<<" "<<c2<<" "<<aa2<<" "<<pp2<<"  "<<c3<<" "<<aa3<<" "<<pp3<<" "<<v1<<" "<<v2<<" "<<chi2<<" "<<erx<<endl;
   cout<<"Error:  "<<erx<<endl;
   cout<<"Final Fluo:  "<<F<<endl;
   clock_t end_time = clock();
   double elapsed_time = static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;
   time1<<count+1<<" "<<elapsed_time<<endl;

    // Print the elapsed time
    std::cout << "Execution time: " << elapsed_time << " seconds" << std::endl;
   if (erx<er) {
    k1=c1;    
    a1=aa1;
    p1=pp1;

    k2=c2;   
    a2=aa2;
    p2=pp2;

    k3=c3;  
    a3=aa3;
    p3=pp3;

    u1=v1;
    u2=v2;
    er=erx;
    cout<<"Mejor ajuste"<<endl; 
    nllm=chi2;
   } 
   count++;
   exp.clear();
   exp.seekg(0, ios::beg);
  } while (count<100);

//Generating the final file

 cout<<"Printing final data"<<endl;
  dtt1=0;
   dtt2=0;
   dtt3=0;
   A1=N;
   A2=N;
   A3=N;
   B1=0;
   B2=0;
   B3=0;
  sq=0;
  for( int i = 0; i < pts; i++  )
   {
    b1=Kinetics1 (dtt1, A1, B1, k1);
    b2=Kinetics2 (dtt2, B1, U1, A2, B2, k2);
    b3=Kinetics2 (dtt3, B2, U2, A3, B3, k3);
    exp>>t>>FF;
    F=Fluo (ff, a1, p1, a2, p2, a3, p3, b1, b2, b3);
    dif=FF-F;
    sq=sq+dif*dif;
    sim<<t<<" "<<FF<<" "<<F<<endl;
   }
   chi2=sq/(sigm*sigm);
   eu=sqrt(sq/pts);
  cout<<"Maximum Likelihood Estimation: "<<endl;
  cout<<"k     "<<"a     "<<"p     "<<endl;
  cout<<k1<<" "<<a1<<" "<<p1<<endl;
  cout<<k2<<" "<<a2<<" "<<p2<<"  "<<u1<<endl;
  cout<<k3<<" "<<a3<<" "<<p3<<" "<<u2<<endl;
  cout<<"chi2: "<<chi2<<endl;
  cout<<"Error: "<<eu<<endl;

 sim.close();
 exp.close();
 par.close();
 time1.close();
 return 0;
}

