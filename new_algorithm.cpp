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

// Total fluorescence
double Fluo (double fo, double a, double b, double c, double d, double e, double f, double f1, double f2, double f3)
{
 double r;
 r=fo+a*f1/(1-b*f1)+c*f2/(1-d*f2)+e*f3/(1-f*f3);
 return r;
}


int main()
{

 clock_t start_time = clock();
 int U1,U2,count,pts, j1, j2;
 double ff,FF, F, dtt,A1,B1,A2,B2,A3,B3,b1,b2,b3,k1,k2,k3,t,tt,a1,a2,a3,p1,p2,p3,dif,sq,erx,u1,u2,v2,chi2, eu;
 double dtt1,dtt2,dtt3;
 double deno1,deno2,deno3, dj_a1, dj_a2, dj_a3, dj_p1, dj_p2, dj_p3, dj_c1, dj_c2, dj_c3, dj_u1, dj_u2, alpha;
 int num1, num2;
 num1 = 1;
 num2 = 5;

 srand(time(NULL));   //required for "randomness"
 ifstream expo ("data_for_analysis.csv");
 ofstream sim;
 sim.open ("model.dat");
 ofstream par ("parmle.dat");
 //ofstream avals ("avals.dat");
 ofstream time1 ("times.dat");
 pts=24000; 

 ff= 1.1882; 

 k1=1.9;    //First kinetic constant (1/ms)
 a1=0.7;
 p1=0.2;

 k2=0.1;   //Second kinetic constant (1/ms)
 a2=0.8;
 p2=0.5;

 k3=0.005;  //Third kinetic constant (1/ms)
 a3=0.7;
 p3=0.5;

 u1=0.9;
 u2=0.8;
 

  par<<"k1     "<<"a1     "<<"p1     "<<"k2     "<<"a2     "<<"p2     "<<"k3     "<<"a3     "<<"p3     "<<"u1   "<<"u2   "<<"nll  "<<endl;

  alpha = 0.09;
  count=0;
  do
  {
   cout<<"Iteration: "<<count+1<<"  ";

   A1=N;
   A2=N;
   A3=N;
   B1=0;
   B2=0;
   B3=0;
   sq=0;
   j1=0;
   j2=0;
   
   
   U1=u1*N;
   U2=u2*N;
   dj_a1 = dj_a2 = dj_a3 = dj_p1 = dj_p2 = dj_p3 = dj_c1 = dj_c2 = dj_c3 = dj_u1 = dj_u2 = 0;

   for( int i = 0; i < pts; i++  )
   {
    
    A1 = int(A1*exp(-0.01*k1));
    B1 = N - A1;
    b1 = B1/N;
    
    if(B1 > U1){
            A2 = int(A2*exp(-0.01*k2));
            j1 = j1 + 1;
            }
    B2 = N - A2;
    b2 = B2/N;
    
    if(B2 > U2){
            A3 = int(A3*exp(-0.01*k3));
            j2 = j2 + 1;
            }
    B3 = N - A3;
    b3 = B3/N;
    
    
    expo>>t>>FF;
    F=Fluo (ff, a1, p1, a2, p2, a3, p3, b1, b2, b3);
    
    //avals<<t<<" "<<A1<<" "<<A2<<" "<<A3<<" "<<F<<endl;
    if((count+1)%num1 == 0){
    sim<<t<<" "<<FF<<" "<<F<<endl;
    }
     
    if((count+1)%num2 == 0){
    sim<<t<<" "<<FF<<" "<<F<<endl;
    }
    dif=F-FF;
    deno1 = 1-(p1*b1);
    deno2 = 1-(p2*b2);
    deno3 = 1-(p3*b3);
    
    dj_a1 = dj_a1 + dif*b1/deno1;
    dj_a2 = dj_a2 + dif*b2/deno2;
    dj_a3 = dj_a3 + dif*b3/deno3;
    
    dj_p1 = dj_p1 + dif*a1*b1*b1/(deno1*deno1);
    dj_p2 = dj_p2 + dif*a2*b2*b2/(deno2*deno2);
    dj_p3 = dj_p3 + dif*a3*b3*b3/(deno3*deno3);
    
    dj_c1 = dj_c1 + dif * a1/(deno1*deno1) * 0.01* A1*exp(-0.01*k1)/N;
    dj_c2 = dj_c2 + dif * a2/(deno2*deno2) * 0.01* A2*exp(-0.01*k2)/N;
    dj_c3 = dj_c3 + dif * a3/(deno3*deno3) * 0.01* A3*exp(-0.01*k3)/N;
    
    if(j1 == 1){
      dj_u1 = dif*a2/(deno2*deno2);
    }
    
    if(j2 == 1){
      dj_u2 = dif*a3/(deno3*deno3);
    }
    
    sq=sq+dif*dif;
   }
   chi2=sq/(sigm*sigm);
   erx=sqrt(sq/pts);
   dj_a1 = dj_a1/pts;
   dj_a2 = dj_a2/pts;
   dj_a3 = dj_a3/pts;
   dj_c1 = dj_c1/pts;
   dj_c2 = dj_c2/pts;
   dj_c3 = dj_c3/pts;
   dj_p1 = dj_p1/pts;
   dj_p2 = dj_p2/pts;
   dj_p3 = dj_p3/pts;
   dj_u1 = dj_u1/pts;
   dj_u2 = dj_u2/pts;
	 
   par<<k1<<" "<<a1<<" "<<p1<<" "<<k2<<" "<<a2<<" "<<p2<<"  "<<k3<<" "<<a3<<" "<<p3<<" "<<u1<<" "<<u2<<" "<<chi2<<" "<<erx<<endl;
   cout<<"Error:  "<<erx<<endl;
   //avals<<endl;
   //avals<<endl;
   //avals<<endl;
   //cout<<"DJ_C1:  "<<dj_c1<<"DJ_C2:  "<<dj_c2<<"DJ_C3:  "<<dj_c3<<endl;
   //
   if((count+1)%num1 == 0){
    num1 = num1*10;
    }

   if((count+1)%num2 == 0){
    num2 = num2*10;
    }

    k1 = k1 - alpha*dj_c1;    
    a1 = a1 - alpha*dj_a1;
    p1 = p1 - alpha*dj_p1;

    k2 = k2- alpha*dj_c2;   
    a2 = a2- alpha*dj_a2;
    p2 = p2- alpha*dj_p2;

    k3 = k3- alpha*dj_c3;  
    a3 = a3- alpha*dj_a3;
    p3 = p3- alpha*dj_p3;
    
    u1 = u1 - alpha*dj_u1;
    u2 = u2 - alpha*dj_u2;

   //clock_t end_time = clock();
   //double elapsed_time = static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;
   //time1<<count+1<<" "<<elapsed_time<<endl;
 
   count++;
   expo.clear();
   expo.seekg(0, ios::beg);
  } while (count<10000);
  
  
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
   U1=u1*N;
   U2=u2*N;
  sq=0;
  j1 = 0;
  j2 = 0;
  for( int i = 0; i < pts; i++  )
   {
    A1 = int(A1*exp(-0.01*k1));
    B1 = N - A1;
    b1 = B1/N;

    if(B1 > U1){
            A2 = int(A2*exp(-0.01*k2));
            j1 = j1 + 1;
            }
    B2 = N - A2;
    b2 = B2/N;

    if(B2 > U2){
            A3 = int(A3*exp(-0.01*k3));
            j2 = j2 + 1;
            }
    B3 = N - A3;
    b3 = B3/N;


    expo>>t>>FF;
    F=Fluo (ff, a1, p1, a2, p2, a3, p3, b1, b2, b3);
    dif=FF-F;
    sq=sq+dif*dif;
    sim<<t<<" "<<FF<<" "<<F<<endl;
    //avals<<t<<" "<<A1<<" "<<A2<<" "<<A3<<" "<<F<<endl;
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
  
  cout<<"djc   "<<"dja    "<<"djp   "<<endl;
  cout<<dj_c1<<" "<<dj_a1<<" "<<dj_p1<<endl;
  cout<<dj_c2<<" "<<dj_a2<<" "<<dj_p2<<endl;
  cout<<dj_c3<<" "<<dj_a2<<" "<<dj_p3<<endl;
  cout<<dj_u1<<" "<<dj_u2<<endl;

 sim.close();
 expo.close();
 par.close();
 time1.close();
 //avals.close();
 return 0;
}

