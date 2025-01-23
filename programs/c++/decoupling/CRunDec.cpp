/*
  CRunDec.cpp, v2.1

  Authors: Florian Herren and Matthias Steinhauser (Jun 2016)

  Changes since v2.0:
  Jun 2016: improved 4-loop coefficients in mMS2mOS, mOS2mSI, mOS2mMS
*/

/*
  Based on:

  CRunDec.cpp, v2.0

  Authors: Florian Herren and Matthias Steinhauser (Jan 2016)

  Changes since v1.1:
  Sep 2015: extended mOS2mMS, mMS2mOS to four loops (see arXiv:1502.01030);  
            mMS2mMS to five loops (see arXiv:1402.6611); 
            DecMqUp/DownOS to four loops (see hep-ph/0512058, hep-ph/0512060
            and arXiv:1502.04719); 
            prepared AsmMSrunexact, AlphasExact for five loops
  Oct 2015: implemented PS, 1S and RS masses, now requires C++11 for tgamma, 
            extended mOS2mSI to four loops (see arXiv:1502.01030)
  Dec 2015: bugfix: nl -> nlq in fas6to5ms, ported DecAsDownSI/DecAsUpSI from RunDec
  Jan 2016: Add RS' mass relations (see arXiv:1502.01030)
*/

/*
  Based on:

  CRunDec.cpp, v1.1

  Author: Barbara Schmidt

  For documentation see

  CRunDec: a C++ package for running and decoupling of the
  strong coupling and quark masses

  by

  Barbara Schmidt, Matthias Steinhauser

  Comput.Phys.Commun. 183 (2012) 1845-1848
  arXiv:1201.6149 
  SFB/CPP-12-03
  TTP12-02

  See also:
  K.G. Chetyrkin, J.H. Kuhn and M. Steinhauser,
  "RunDec: A Mathematica package for running and decoupling of the strong
  coupling and quark masses"
  Comput. Phys. Commun.  133 (2000) 43,
  arXiv:hep-ph/0004189

  May 2013: minor modifications to avoid "-Wunused-parameter"
            during compilation with g++ -Wall -Wextra
  Sep 2013: bug fix in fRungeKuttaImpl (thanks to Stephen Jones)
  Oct 2014: bug fix in AsmMSrunexact
*/

/*
License:

CRunDec is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include "CRunDec.h"

#define BETA4 100.0

using namespace std;

// Default constructor:
CRunDec::CRunDec(){
     for(int i=0; i<4; i++){
       mq[i]=0.;
       nfMmu[i].Mth=0.;
       nfMmu[i].muth=0.;
       nfMmu[i].nf=0;
     }
     AM.Asexact=0.;
     AM.mMSexact=0.;
}

// Constructor called with number for active flavours nf
CRunDec::CRunDec(int n){
     double nf=(double) n;
     Nf=nf;
     Beta[0]= (double)0.25*(11.0 - 2.*nf/3.0);
     Beta[1]= (double)(102.0 - 38.*nf/3.0)/16.0;
     Beta[2]= (double)(0.5*2857. - 5033.*nf/18.0 + 325.*nf*nf/54.0)/64.0;
     Beta[3]= (double)(149753./6.0 + 3564.*Zeta3 + 
              (-1078361./162.0 - 6508.*Zeta3/27.0)*nf + 
              (50065./162.0 + 6472.*Zeta3/81.0)*nf*nf + 
              1093.*nf*nf*nf/729.0)/256.0;
     Beta[4]= (double)0.25*(11.0 - 2.*nf/3.0)*BETA4;

     Gamma[0]=(double)1.0;
     Gamma[1]=(double)(202./3.-20.*nf/9.)/16.;  
     Gamma[2]=(double)(1249. + (-2216./27. - 160.*Zeta3/3.)*nf-
              140.*nf*nf/81.)/64.;
     Gamma[3]=(double)(4603055./162. + 135680.*Zeta3/27. - 8800.*Zeta5 +
              (-91723./27. - 34192.*Zeta3/9. + 
              880.*Zeta4 + 18400.*Zeta5/9.)*nf +
              (5242./243. + 800.*Zeta3/9. - 160.*Zeta4/3.)*nf*nf +
              (-332./243. + 64.*Zeta3/27.)*nf*nf*nf)/256.;
     Gamma[4]=(double)(99512327./162. + 46402466.*Zeta3/243. + 96800.*Zeta3*Zeta3 - 698126.*Zeta4/9.
              -231757160.*Zeta5/243. + 242000.*Zeta6 + 412720.*Zeta7
              +nf*(-150736283./1458. - 12538016.*Zeta3/81. - 75680.*Zeta3*Zeta3/9. + 2038742.*Zeta4/27.
              + 49876180.*Zeta5/243. - 638000.*Zeta6/9. - 1820000.*Zeta7/27.)
              +nf*nf*(1320742./729. + 2010824.*Zeta3/243. + 46400.*Zeta3*Zeta3/27. - 166300.*Zeta4/27. - 264040.*Zeta5/81. + 92000.*Zeta6/27.)
              +nf*nf*nf*(91865./1458. + 12848.*Zeta3/81. + 448.*Zeta4/9. - 5120.*Zeta5/27.)
              +nf*nf*nf*nf*(-260./243. - 320.*Zeta3/243. + 64.*Zeta4/27.)
              )/(4*4*4*4*4);
     

     for(int i=0; i<5; i++) {
       B[i]=Beta[i]/Beta[0];
       C[i]=Gamma[i]/Beta[0];
     }
     for(int i=0; i<4; i++) {
       mq[i]=0.;
       nfMmu[i].Mth=0.;
       nfMmu[i].muth=0.;
       nfMmu[i].nf=0;
     } 
     AM.Asexact=0.;
     AM.mMSexact=0.;        
}

// Define constants (if not already done in constructor)
void CRunDec::SetConstants(int n){
     double nf=(double) n;
     if(Nf!=nf){
     }
     Nf=nf;
     Beta[0]= (double)0.25*(11.0 - 2.*nf/3.0);
     Beta[1]= (double)(102.0 - 38.*nf/3.0)/16.0;
     Beta[2]= (double)(0.5*2857. - 5033.*nf/18.0 + 325.*nf*nf/54.0)/64.0;
     Beta[3]= (double)(149753./6.0 + 3564.*Zeta3 + 
              (-1078361./162.0 - 6508.*Zeta3/27.0)*nf + 
              (50065./162.0 + 6472.*Zeta3/81.0)*nf*nf + 
              1093.*nf*nf*nf/729.0)/256.0;
     Beta[4]= (double)0.25*(11.0 - 2.*nf/3.0)*BETA4;
     
     Gamma[0]=(double)1.0;
     Gamma[1]=(double)(202./3.-20.*nf/9.)/16.; 
     Gamma[2]=(double)(1249. + (-2216./27. - 160.*Zeta3/3.)*nf-
              140.*nf*nf/81.)/64.;
     Gamma[3]=(double)(4603055./162. + 135680.*Zeta3/27. - 8800.*Zeta5 +
              (-91723./27. - 34192.*Zeta3/9. + 
              880.*Zeta4 + 18400.*Zeta5/9.)*nf +
              (5242./243. + 800.*Zeta3/9. - 160.*Zeta4/3.)*nf*nf +
              (-332./243. + 64.*Zeta3/27.)*nf*nf*nf)/256.;
     Gamma[4]=(double)(99512327./162. + 46402466.*Zeta3/243. + 96800.*Zeta3*Zeta3 - 698126.*Zeta4/9.
              -231757160.*Zeta5/243. + 242000.*Zeta6 + 412720.*Zeta7
              +nf*(-150736283./1458. - 12538016.*Zeta3/81. - 75680.*Zeta3*Zeta3/9. + 2038742.*Zeta4/27.
              + 49876180.*Zeta5/243. - 638000.*Zeta6/9. - 1820000.*Zeta7/27.)
              +nf*nf*(1320742./729. + 2010824.*Zeta3/243. + 46400.*Zeta3*Zeta3/27. - 166300.*Zeta4/27. - 264040.*Zeta5/81. + 92000.*Zeta6/27.)
              +nf*nf*nf*(91865./1458. + 12848.*Zeta3/81. + 448.*Zeta4/9. - 5120.*Zeta5/27.)
              +nf*nf*nf*nf*(-260./243 - 320.*Zeta3/243. + 64.*Zeta4/27.)
              )/(4*4*4*4*4);

       
     for(int i=0; i<5; i++) {
       B[i]=Beta[i]/Beta[0];
       C[i]=Gamma[i]/Beta[0];
     }                          
}


// Function int CRunDec::GetNf()
// Returns number of active flavours.
int CRunDec::GetNf(){
     return (int)Nf;
}

// Function void CRunDec::SetNf(int nf)
// Set the private component Nf to the number of active flavours.
void CRunDec::SetNf(int nf){
     this->SetConstants(nf);
}

// Aux. function to exit function.
int CRunDec::Abbruch(void){
     RETURN
}

// Function double CRunDec::LamExpl(double AlphaS, double Mu, int nl)
// Compute \Lambda using eq.(4) of [RunDec].
double CRunDec::LamExpl(double AlphaS, double Mu, int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN
     }
     double A=AlphaS/Pi;
     double sum[4];
     sum[0]= 1./(A*Beta[0]);
     sum[1]= (B[1]*log(A))/Beta[0]
              + (B[1]/Beta[0])*log(Beta[0]);
     sum[2]= (B[2]*A - B[1]*B[1]*A )/Beta[0];
     sum[3]= (0.5*B[3]*A*A-B[1]*B[2]*A*A+ 0.5*B[1]*B[1]*B[1]*A*A)/Beta[0];

     double LogM2L2=0.0;
     for(int i=1; i<=nl; i++){
       LogM2L2+=sum[i-1];
     }
               
     double Lambda= Mu*exp(-0.5*LogM2L2);
     return Lambda;
}

// Function double CRunDec::AlphasLam(double Lambda, double Mu, int nl)
// Compute \alpha_s using eq.(5) of [RunDec]. 
double CRunDec::AlphasLam(double Lambda, double Mu, int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(Mu/Lambda<1.5){
       cout<<"WARNING: the ratio \\mu/\\lambda = "<< Mu/Lambda
             <<" is very small!"<<endl; 
       RETURN  
     }
     double L=log(Mu*Mu/(Lambda*Lambda));
     double h=1/(L*Beta[0]);
     double c=log(L);
     double sum[4];
     sum[0]= h;
     sum[1]= -h*h*B[1]*c;
     sum[2]= + h*h*h*(B[1]*B[1]*(c*c-c-1)+B[2]);
     sum[3]= h*h*h*h*(B[1]*B[1]*B[1]*(-c*c*c+2.5*c*c+2*c-0.5) 
                   -3*B[1]*B[2]*c +0.5*B[3]);
     double a=0.0;
     for(int i=1; i<=nl; i++){
       a+=sum[i-1];
     }
     return a*Pi;
}

// Eq.(5) rewritten in a form suitable to determine zero.
double CRunDec::fSetAsL(double Lambda, double Mu, int nl, double AlphaS){
     double L=log(Mu*Mu/(Lambda*Lambda));
     double h=1/(L*Beta[0]);
     double c=log(L);
     double sum[4];
     sum[0]= h;
     sum[1]= -h*h*B[1]*c;
     sum[2]= + h*h*h*(B[1]*B[1]*(c*c-c-1)+B[2]);
     sum[3]= h*h*h*h*(B[1]*B[1]*B[1]*(-c*c*c+2.5*c*c+2*c-0.5) 
                   -3*B[1]*B[2]*c +0.5*B[3]);
     double Add=0.0;
     for(int i=1; i<=nl; i++) Add+=sum[i-1];
     return (Add-(AlphaS/Pi));
}

// Function double CRunDec::LamImpl(double AlphaS, double Mu,int nl)
// Compute \Lambda using eq.(5) of [RunDec]. 
double CRunDec::LamImpl(double AlphaS, double Mu,int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double epsilonX= 1e-8;
     double Lambda0=LamExpl(AlphaS,Mu,nl);
     double x0=Lambda0 - 0.2*Lambda0;
     double x1=Lambda0 + 0.2*Lambda0;
     double f0= this->fSetAsL(x0,Mu,nl,AlphaS);
     double f1= this->fSetAsL(x1,Mu,nl,AlphaS);
     if(f0*f1>0){
       cout<<"WARNING: No root can be calculatet!"<<endl;
       RETURN
     }
     double xTest;
     double fTest;
     do{
       xTest= (x0+x1)/2;
       fTest= fSetAsL(xTest,Mu,nl,AlphaS);
       if(f0*fTest<0){x1= xTest;}
       else {x0= xTest;}
     }
     while(abs(x1-x0)>= epsilonX);
     double Lambda=xTest;
     return Lambda;
}

// Right-hand side of differential equation for \alpha_s times 1/\mu^2,
// see eq.(1) of [1].
double fSetdydx(CRunDec S, double A,int nl){ 
     double f=0.0;
     double sum[5];
     double B=A*A;
     sum[0]=-S.Beta[0]*B;    
     sum[1]=-S.Beta[1]*B*A;
     sum[2]=-S.Beta[2]*B*B;
     sum[3]=-S.Beta[3]*B*B*A; 
     sum[4]=-S.Beta[4]*B*B*B;   
     for(int i=1; i<= nl; i++) {
       f+=sum[i-1];
     }
     return (f*2);
}

// Implicit Runge-Kutte step (4th order)
// Call with x=\mu^2, y=\alpha_s(\mu)/\pi, step size h, number of loops nl 
// New y value is returned at x+h
double CRunDec::fRungeKuttaImpl(double &x, double y,double &htry, int nl, 
                              double (*f)(CRunDec,double, int)){
     // Precision
     double eps=1e-10;
     double yerr,ytemp,htemp, hnext;
     double h=htry;
     double xnew; // new variable   
     double k1,k2,k3,k4,k5,k6;
     for(;;){
       k1=h*f(*this,y,nl);
       k2=h*f(*this,y+b21*k1,nl);
       k3=h*f(*this,y+b31*k1+b32*k2,nl);
       k4=h*f(*this,y+b41*k1+b42*k2+b43*k3,nl);
       k5=h*f(*this,y+b51*k1+b52*k2+b53*k3+b54*k4,nl);
       k6=h*f(*this,y+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5,nl);
       // y value at x+h as a sum of the previous value and the
       // correspondingly weighted function evaluations
       ytemp= y+ c1*k1+ c2*k2+ c3*k3+ c4*k4+ c5*k5+ c6*k6;
       // Estimate of uncertainty
       yerr=dc1*k1 + dc2*k2 + dc3*k3 + dc4*k4 + dc5*k5 + dc6*k6;
       double err=0.;
       err=fmax(err,fabs(yerr/eps));
       
       // Uncertainty too big? -> Discard result and reduce step size
       if(err>1.){      
         htemp=0.9*h*pow(err,-0.25);
         if(h>=0.){h=fmax(htemp,0.1*h);}
         else{h=fmin(htemp,0.1*h);}
         xnew=x+h;  // modification to previous code
         //decide whether reduced stepsize is still big enough 
         //(in order to prevent a closed loop)
         if(xnew==x){cout<<"stepsize too small"<<endl; RETURN} 
         continue;          
       }
       else{
         if(err>1.89e-4){
         hnext=0.9*h*pow(err,-0.2);
         }
         // Uncertainty OK? -> take y value, increase h
         else{
           hnext=5.*h;
         }
         x+=h;
         
         y=ytemp;
         htry=hnext;
         break;
       }       
     }
     return y; 
}

// Function: double CRunDec::AlphasExact(double AlphaS0, double Mu0, 
//                           double MuEnd, int nl)
// Compute \alpha_s using eq.(1) of [RunDec]
double CRunDec::AlphasExact(double AlphaS0, double Mu0, double MuEnd, int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN
     }
     double Lambda=LamExpl(AlphaS0,Mu0,(nl < 5) ? nl : 4);
     if(MuEnd/Lambda<1.5){
       cout<<"WARNING: the ratio \\mu/\\lambda = "<< MuEnd/Lambda
         <<" is very small!"<<endl;
       RETURN
     }
     double x,y;
     x=log(Mu0);
     y=AlphaS0/Pi;
     double h;
     
     if(Mu0<MuEnd){
       h=1e-4;
       while(x<(log(MuEnd))){
         y=this->fRungeKuttaImpl(x,y,h,nl,fSetdydx);
         if(x+h>=log(MuEnd)){
           h=log(MuEnd)-x;
         }
       }
     return(y*Pi);
     }
     else{h=-1e-4;}
     while(x>(log(MuEnd))){
       y=this->fRungeKuttaImpl(x,y,h,nl,fSetdydx);
       if(x+h<=log(MuEnd)){
         h=log(MuEnd)-x;
       }
     }
     return(y*Pi);
} 

// Eq.(10) of [RunDec]
double CRunDec::fSetcx(double x, int nl){
     double sum[5];
     sum[0]=1;
     sum[1]=(C[1]-B[1]*C[0])*x;
     sum[2]=0.5*((C[1]-B[1]*C[0])*(C[1]-B[1]*C[0]) + C[2] - B[1]*C[1]+
           B[1]*B[1]*C[0] - B[2]*C[0])*x*x;
     sum[3]=((C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])/6. +
           0.5*(C[1]-B[1]*C[0])*(C[2]-B[1]*C[1]+B[1]*B[1]*C[0]-B[2]*C[0])+
           (C[3]-B[1]*C[2]+B[1]*B[1]*C[1]-B[2]*C[1]-B[1]*B[1]*B[1]*C[0] +
            2.*B[1]*B[2]*C[0] - B[3]*C[0])/3.)*x*x*x;
     sum[4]=((C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])/24.
             +(C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])*(C[2]/2. - B[1]*C[1]/2.+ B[1]*B[1]*C[0]/2. - B[2]*C[0]/2.)/2.
             +(C[2]/2. - B[1]*C[1]/2.+ B[1]*B[1]*C[0]/2. - B[2]*C[0]/2.)*(C[2]/2. - B[1]*C[1]/2.+ B[1]*B[1]*C[0]/2. - B[2]*C[0]/2.)/2.
             +((C[1]-B[1]*C[0])*(C[3]-B[1]*C[2]+B[1]*B[1]*C[1]-B[2]*C[1]-B[1]*B[1]*B[1]*C[0] + 2.*B[1]*B[2]*C[0] - B[3]*C[0])/3.)
             +B[1]*B[1]*B[1]*B[1]*C[0]/4. - 3.*B[1]*B[1]*B[2]*C[0]/4. + B[2]*B[2]*C[0]/4. + B[1]*B[3]*C[0]/2. - B[4]*C[0]/4.
             -B[1]*B[1]*B[1]*C[1]/4. + B[1]*B[2]*C[1]/2. - B[3]*C[1]/4. + B[1]*B[1]*C[2]/4. - B[2]*C[2]/4. - B[1]*C[3]/4. + C[4]/4.)*x*x*x*x;
     double erg=0.0;
     for(int i=1; i<=nl; i++){
       erg+=sum[i-1];
     }
     return (pow(x,C[0])*erg);             
}

// Function double CRunDec::mMS2mMS(double Mu0, double AlphaS0, 
//                          double AlphaSEnd, int nl)
// Compute m_q(\mu) using eqs.(9) and (10) of [RunDec]
double CRunDec::mMS2mMS(double Mu0, double AlphaS0, double AlphaSEnd, int nl){
     if(nl<0||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl==0){
       return Mu0;
     }
     double cAlphaS0= this->fSetcx(AlphaS0/Pi, nl);
     double cAlphaSEnd= this->fSetcx(AlphaSEnd/Pi, nl);
     return Mu0*cAlphaSEnd/cAlphaS0;        
}

// Aux. functions (r.h.s of diff. eqs. for alpha_s and m_q)
double fSetdydxM1(CRunDec S,double A, double M){
     return (M*(S.Gamma[0])/(S.Beta[0]*A));
}

double fSetdydxa1(CRunDec S,double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A));      
}

double fSetdydxM2(CRunDec S,double A, double M){
     return (M*(S.Gamma[0]+S.Gamma[1]*A)/(S.Beta[0]*A+S.Beta[1]*A*A));   
}

double fSetdydxa2(CRunDec S,double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A+S.Beta[1]*A*A*A));      
}

double fSetdydxM3(CRunDec S,double A, double M){
     return (M*(S.Gamma[0]+S.Gamma[1]*A+S.Gamma[2]*A*A)/
           (S.Beta[0]*A+S.Beta[1]*A*A+S.Beta[2]*A*A*A));     
}

double fSetdydxa3(CRunDec S,double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A+S.Beta[1]*A*A*A+S.Beta[2]*A*A*A*A));      
}

double fSetdydxM4(CRunDec S,double A, double M){
     return (M*(S.Gamma[0]+S.Gamma[1]*A+S.Gamma[2]*A*A+S.Gamma[3]*A*A*A)/
           (S.Beta[0]*A+S.Beta[1]*A*A+S.Beta[2]*A*A*A+S.Beta[3]*A*A*A*A)); 
}

double fSetdydxa4(CRunDec S,double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A+S.Beta[1]*A*A*A+S.Beta[2]*A*A*A*A+
                S.Beta[3]*A*A*A*A*A));      
}

double fSetdydxM5(CRunDec S,double A, double M){
     return (M*(S.Gamma[0]+S.Gamma[1]*A+S.Gamma[2]*A*A+S.Gamma[3]*A*A*A + S.Gamma[4]*A*A*A*A)/
           (S.Beta[0]*A+S.Beta[1]*A*A+S.Beta[2]*A*A*A+S.Beta[3]*A*A*A*A+S.Beta[4]*A*A*A*A*A)); 
}

double fSetdydxa5(CRunDec S,double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A+S.Beta[1]*A*A*A+S.Beta[2]*A*A*A*A+
                S.Beta[3]*A*A*A*A*A+S.Beta[4]*A*A*A*A*A*A));      
}  

// Runge-Kutta step for implicit procedure
double CRunDec::fRKSchritt(double x,double y,double h,double &yerr,
                    double (*f)(CRunDec, double ,double)){
     double k1,k2,k3,k4,k5,k6; 
     k1=h*f(*this,x,y);
     k2=h*f(*this,x+a2*h,y+b21*k1);
     k3=h*f(*this,x+a3*h,y+b31*k1+b32*k2);
     k4=h*f(*this,x+a4*h,y+b41*k1+b42*k2+b43*k3);
     k5=h*f(*this,x+a5*h,y+b51*k1+b52*k2+b53*k3+b54*k4);
     k6=h*f(*this,x+a6*h,y+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5);

     yerr=(dc1*k1 + dc2*k2 + dc3*k3 + dc4*k4 + dc5*k5 + dc6*k6);
     return (y+ c1*k1+ c2*k2+ c3*k3+ c4*k4+ c5*k5+ c6*k6);      
}  

// Function AsmMS CRunDec::AsmMSrunexact(double mMu, double AlphaS0, double Mu0, 
//                         double MuEnd, int nl)
// Compute \alpha_s and m_q solving diff. eqs. simultaneously.
AsmMS CRunDec::AsmMSrunexact(double mMu, double AlphaS0, double Mu0, 
                             double MuEnd, int nl){
     AsmMS Erg;
     if(nl<0||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       this->Abbruch();
     }
     if(nl==0){
       Erg.Asexact=AlphaS0;
       Erg.mMSexact=mMu;
       return Erg;
     }
     double yerr0,ytemp0,xnew,h=1e-3;
     float eps=1e-15;
     float errmax;
     double x0=log(Mu0);
     double y0=AlphaS0/Pi;
     double xEnd=log(MuEnd);
     double yscal0=abs(x0)+abs(h*y0);
     
     double (*falpha)(CRunDec, double ,double);
     double (*fmMS)(CRunDec, double ,double);
     
     if(nl==1){
       falpha=fSetdydxa1;
       fmMS=fSetdydxM1;
     }
     if(nl==2){
       falpha=fSetdydxa2;
       fmMS=fSetdydxM2;
     } 
     if(nl==3){
       falpha=fSetdydxa3;
       fmMS=fSetdydxM3;
     }
     if(nl==4){
       falpha=fSetdydxa4;
       fmMS=fSetdydxM4;
     }
     if(nl==5){
       falpha=fSetdydxa5;
       fmMS=fSetdydxM5;
     }
     
     if(Mu0<MuEnd){
       h=1e-2;
       while(x0<xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,falpha);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;  
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h>=xEnd){
           h=xEnd-x0;
         }//if
       }//while
       Erg.Asexact=y0*Pi;   
       
       x0=AlphaS0/Pi;
       xEnd=y0;
       y0=mMu;
       yscal0=abs(x0)+abs(h*y0);
       eps=1e-10;
       h=-1e-3;

       while(x0>xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,fmMS);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h<=xEnd){
           h=xEnd-x0;
         }//if
       }//while  
       Erg.mMSexact=y0;
       return Erg;
     }//if
     else{
       h=-1e-2;
       while(x0>xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,falpha);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h<=xEnd){
           h=xEnd-x0;
         }//if
       }//while
       Erg.Asexact=y0*Pi; 

       x0=AlphaS0/Pi;     
       xEnd=y0;
       y0=mMu;
       yscal0=abs(x0)+abs(h*y0);
       eps=1e-10;
       h=1e-3;
       while(x0<xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,fmMS);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h>=xEnd){
           h=xEnd-x0;
         }//if
       }//while
       Erg.mMSexact=y0;
       return Erg;
     }//else
}   

// Coefficients of eq.(13) of [RunDec]
double CRunDec::fMsFromOs1(double mu, double M){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) (-cf - (3.*cf*lmM)/4.);
     return erg;
}

double CRunDec::fMsFromOs2(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) ((-1111.*ca*cf)/384. + (7.*cf*cf)/128. -
     (185.*ca*cf*lmM)/96. + (21.*cf*cf*lmM)/32. - (11.*ca*cf*lmM*lmM)/32. +
     (9.*cf*cf*lmM*lmM)/32. + (143.*cf*tr)/96. + (13.*cf*lmM*tr)/24. + 
     (cf*lmM*lmM*tr)/8. +
     (71.*cf*nl*tr)/96. + (13.*cf*lmM*nl*tr)/24. + (cf*lmM*lmM*nl*tr)/8. +
     (ca*cf*Zeta2)/2 - (15.*cf*cf*Zeta2)/8. - (3.*ca*cf*log(2)*Zeta2)/2. + 
     3.*cf*cf*log(2)*Zeta2 -
     cf*tr*Zeta2 + (cf*nl*tr*Zeta2)/2. + (3.*ca*cf*Zeta3)/8. - 
     (3.*cf*cf*Zeta3)/4.);
     return erg;
}

double CRunDec::fMsFromOs3(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) (lmM*lmM*(-2341.*ca*ca*cf + 1962.*ca*cf*cf - 243.*cf*cf*cf 
     + 1492.*ca*cf*tr -
      468.*cf*cf*tr + 1492.*ca*cf*nl*tr - 468.*cf*cf*nl*tr - 208.*cf*tr*tr -
      416.*cf*nl*tr*tr - 208.*cf*nl*nl*tr*tr))/1152. +
     (lmM*lmM*lmM*(-242.*ca*ca*cf + 297.*ca*cf*cf - 81.*cf*cf*cf + 
     176.*ca*cf*tr - 108.*cf*cf*tr + 176.*ca*cf*nl*tr - 108.*cf*cf*nl*tr - 
     32.*cf*tr*tr - 64.*cf*nl*tr*tr - 32.*cf*nl*nl*tr*tr))/1152. +
     (lmM*(-105944.*ca*ca*cf + 52317.*ca*cf*cf - 13203.*cf*cf*cf + 
     74624.*ca*cf*tr -
     5436.*cf*cf*tr + 55616.*ca*cf*nl*tr + 2340.*cf*cf*nl*tr - 
     12608.*cf*tr*tr -
     18304.*cf*nl*tr*tr - 5696.*cf*nl*nl*tr*tr + 12672.*ca*ca*cf*Zeta2 -
     52704.*ca*cf*cf*Zeta2 + 19440.*cf*cf*cf*Zeta2 - 
     38016.*ca*ca*cf*log(2)*Zeta2 +
     91584.*ca*cf*cf*log(2)*Zeta2 - 31104.*cf*cf*cf*log(2)*Zeta2 - 
     29952.*ca*cf*tr*Zeta2 +
     27648.*cf*cf*tr*Zeta2 + 13824.*ca*cf*log(2)*tr*Zeta2 - 
     27648.*cf*cf*log(2)*tr*Zeta2 +
     8064.*ca*cf*nl*tr*Zeta2 + 12096.*cf*cf*nl*tr*Zeta2 + 
     13824.*ca*cf*log(2)*nl*tr*Zeta2 -
     27648.*cf*cf*log(2)*nl*tr*Zeta2 + 9216.*cf*tr*tr*Zeta2 + 
     4608.*cf*nl*tr*tr*Zeta2 -
     4608.*cf*nl*nl*tr*tr*Zeta2 + 9504.*ca*ca*cf*Zeta3 - 
     22896.*ca*cf*cf*Zeta3 +
     7776.*cf*cf*cf*Zeta3 + 6912.*ca*cf*tr*Zeta3 - 3456.*cf*cf*tr*Zeta3 +
     6912.*ca*cf*nl*tr*Zeta3 - 3456.*cf*cf*nl*tr*Zeta3))/13824.;
     return erg;
}

//z[3,m](M) according to eq.(15)
double CRunDec::fZmM(double nl){
     double erg;
     erg= (double) -9478333./93312. + 55.*log(2)*log(2)*log(2)*log(2)/162. +
            (-644201./6480. + 587.*log(2)/27. + 44.*log(2)*log(2)/27.)*Zeta2 -
            61.*Zeta3/27. + 3475*Zeta4/432. + 1439.*Zeta2*Zeta3/72. -
            1975.*Zeta5/216. + 220.*A4/27. + nl*(246643./23328. - 
            log(2)*log(2)*log(2)*log(2)/81. +(967./108. + 22.*log(2)/27. -
            4.*log(2)*log(2)/27.)*Zeta2 + 241.*Zeta3/72. - 305.*Zeta4/108. -
            8.*A4/27.) + nl*nl*(-2353./23328. - 13.*Zeta2/54 - 7.*Zeta3/54.);
     return erg;
}

double CRunDec::fMsFromOs4(double mu, double M, double nf, double err){
     double lmM=log((mu*mu)/(M*M));
     double erg;

     erg =  - 3654.15040757339*err - 1524.2292266911543*lmM - 288.778291935394*lmM*lmM - 32.54735725308642*lmM*lmM*lmM - 
            1.85546875*lmM*lmM*lmM*lmM +
            (-1. + nf)*(-1. + nf)*(-1. + nf)*(0. + 0.678141025604516*err + 0.3205521521864135*lmM + 0.0800290327210927*lmM*lmM +      
            0.010030864197530864*lmM*lmM*lmM) + 
            (-1. + nf)*(-1. + nf)*(0. - 43.48241924867489*err - 19.82672048099557*lmM - 4.482957520194182*lmM*lmM - 
            0.5270061728395061*lmM*lmM*lmM - 0.04108796296296297*lmM*lmM*lmM*lmM) +
            (-1. + nf)*(0. + 756.9421565599532*err + 330.1770776731065*lmM + 67.99849534415492*lmM*lmM + 
            7.595293209876542*lmM*lmM*lmM + 0.48119212962962954*lmM*lmM*lmM*lmM);
     return erg;
}

// Function: double CRunDec::mOS2mMS(double mOS, double mq[], double asmu,
//                           double Mu,int nl)
double CRunDec::mOS2mMS(double mOS, double mq[], double asmu,double Mu,int nl,double err){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN 
     }
     if(nl == 4 && (Nf<4||Nf>6)){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR 4 LOOPS and "<< nl << " FLAVORS"<<endl;
       RETURN  
     }
     double sum[5];
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fMsFromOs1(Mu, mOS))/Pi;
     sum[2]=asmu*asmu*((this-> fMsFromOs2(Mu, mOS, Nf-1))
           -4.*(this->fDelta(mOS,mq)/3.))/(Pi*Pi); 
     sum[3]=asmu*asmu*asmu*(this-> fMsFromOs3(Mu, mOS,Nf-1)+
     this->fZmM(Nf-1))/(Pi*Pi*Pi);
     sum[4] = asmu*asmu*asmu*asmu*(this->fMsFromOs4(Mu,mOS,Nf,err))/(Pi*Pi*Pi*Pi); ;
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     return mOS*erg;       
}

// Function: double CRunDec::mMS2mSI(double mMS, double asmu, double mu, int nl)
double CRunDec::mMS2mSI(double mMS, double asmu, double mu, int nl){   
  double epsilonX = 1e-8;
  AsmMS asmq;
  asmq.Asexact  = asmu;
  asmq.mMSexact = mMS;
  for (;;) {
    double mbold = asmq.mMSexact;
    asmq = AsmMSrunexact(mMS, asmu, mu, mbold, nl);
    if (abs(asmq.mMSexact - mbold) < epsilonX) break;
  }
  return asmq.mMSexact;
}

// Coefficients of eq.(18) of [RunDec]
double CRunDec::fMsFromRi1(void){
     return (double) -4./3.;
}
double CRunDec::fMsFromRi2(void){
     double erg= (double) -995./72. + 19.*Zeta3/6. + 89.*Nf/144.;
     return erg;       
}
double CRunDec::fMsFromRi3(void){
     double erg= (double) -6663911./41472. + 408007.*Zeta3/6912. -185.*Zeta5/36.
                + (118325./7776. + 5*Zeta4/12. - 617.*Zeta3/216.)*Nf +
                (-4459./23328. - Zeta3/54.)*Nf*Nf;
     return erg;                    
}

// Function: double CRunDec::mRI2mMS(double mRI, double asmu, int nl)
double CRunDec::mRI2mMS(double mRI, double asmu, int nl){
     if(nl<0||nl>3){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double sum[4];
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fMsFromRi1())/Pi;
     sum[2]=asmu*asmu*(this-> fMsFromRi2())/(Pi*Pi);
     sum[3]=asmu*asmu*asmu*(this-> fMsFromRi3())/(Pi*Pi*Pi);
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{
       for(int i=0; i<=nl; i++){
        erg+=sum[i];
       }
     }
     return mRI*erg;         
}

// Function: double CRunDec::mMS2mRGI(double mMS, double asmu, int nl)
double CRunDec::mMS2mRGI(double mMS, double asmu, int nl){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl==0){
       return (double) mMS;
     }
     else{    
       double cAsmu= this->fSetcx(asmu/Pi, nl); 
       return (double) mMS/cAsmu;
     }
}

// Function: double CRunDec::mRGI2mMS(double mRGI, double asmu, int nl)
double CRunDec::mRGI2mMS(double mRGI, double asmu, int nl){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN   
     }  
     if(nl==0){
       return (double) mRGI;
     }
     double cAsmu= this->fSetcx(asmu/Pi, nl); 
     return (double) mRGI*cAsmu;             
}

// Coefficients of eq.(17) of [RunDec]
double CRunDec::fOsFromMs1(double mu, double M){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) (cf + (3.*cf*lmM)/4.);
     return erg;
}

double CRunDec::fOsFromMs2(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) ((1111.*ca*cf)/384. - (71.*cf*cf)/128. -
     (143.*cf*tr)/96. - (71.*cf*nl*tr)/96. +
     lmM*((185.*ca*cf)/96. - (9.*cf*cf)/32. - (13.*cf*tr)/24. - 
     (13.*cf*nl*tr)/24.) +
     lmM*lmM*((11.*ca*cf)/32. + (9.*cf*cf)/32. - (cf*tr)/8. - (cf*nl*tr)/8.) -
     (ca*cf*Zeta2)/2. + (15.*cf*cf*Zeta2)/8. + (3.*ca*cf*log(2)*Zeta2)/2. 
     - 3.*cf*cf*log(2)*Zeta2 +
     cf*tr*Zeta2 - (cf*nl*tr*Zeta2)/2. - (3.*ca*cf*Zeta3)/8. + 
     (3.*cf*cf*Zeta3)/4.);
     return erg;
}

double CRunDec::fOsFromMs3(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;

     erg= (double) (lmM*lmM*lmM*((121.*ca*ca*cf)/576. + (33.*ca*cf*cf)/128. +
     (9.*cf*cf*cf)/128. - (11.*ca*cf*tr)/72. - (3.*cf*cf*tr)/32. - 
     (11.*ca*cf*nl*tr)/72. -
     (3.*cf*cf*nl*tr)/32. + (cf*tr*tr)/36. + (cf*nl*tr*tr)/18. +
     (cf*nl*nl*tr*tr)/36.) + lmM*lmM*((2341.*ca*ca*cf)/1152. 
     + (21.*ca*cf*cf)/64. -
     (63.*cf*cf*cf)/128. - (373.*ca*cf*tr)/288. - (3.*cf*cf*tr)/32. -
     (373.*ca*cf*nl*tr)/288. - (3.*cf*cf*nl*tr)/32. + (13.*cf*tr*tr)/72. +
     (13.*cf*nl*tr*tr)/36. + (13.*cf*nl*nl*tr*tr)/72.) +
     lmM*((13243.*ca*ca*cf)/1728. - (4219.*ca*cf*cf)/1536. + 
     (495.*cf*cf*cf)/512. -
     (583.*ca*cf*tr)/108. - (307.*cf*cf*tr)/384. - (869.*ca*cf*nl*tr)/216. -
     (91.*cf*cf*nl*tr)/384. + (197.*cf*tr*tr)/216. + (143.*cf*nl*tr*tr)/108. +
     (89.*cf*nl*nl*tr*tr)/216. - (11.*ca*ca*cf*Zeta2)/12. + 
     (49.*ca*cf*cf*Zeta2)/16. +
     (45.*cf*cf*cf*Zeta2)/32. + (11.*ca*ca*cf*log(2)*Zeta2)/4. - 
     (35.*ca*cf*cf*log(2)*Zeta2)/8. -
     (9.*cf*cf*cf*log(2)*Zeta2)/4. + (13.*ca*cf*tr*Zeta2)/6. - 
     (cf*cf*tr*Zeta2)/2. -
     ca*cf*log(2)*tr*Zeta2 + 2.*cf*cf*log(2)*tr*Zeta2 - 
     (7.*ca*cf*nl*tr*Zeta2)/12. -
     (13.*cf*cf*nl*tr*Zeta2)/8. - ca*cf*log(2)*nl*tr*Zeta2 + 
     2.*cf*cf*log(2)*nl*tr*Zeta2 -
     (2.*cf*tr*tr*Zeta2)/3. - (cf*nl*tr*tr*Zeta2)/3. + 
     (cf*nl*nl*tr*tr*Zeta2)/3. -
     (11.*ca*ca*cf*Zeta3)/16. + (35.*ca*cf*cf*Zeta3)/32. + 
     (9.*cf*cf*cf*Zeta3)/16. -
     (ca*cf*tr*Zeta3)/2. + (cf*cf*tr*Zeta3)/4. - (ca*cf*nl*tr*Zeta3)/2. +
     (cf*cf*nl*tr*Zeta3)/4.));
     return erg;
}

double CRunDec::fOsFromMs4(double mu, double M, double nf, double err){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg = 3567.602784989066*err + 1727.2260148986106*lmM + 409.2429990574718*lmM*lmM + 66.93663194444443*lmM*lmM*lmM + 
           8.056278935185185*lmM*lmM*lmM*lmM + 
           (-1. + nf)*(-1. + nf)*(-1. + nf)*(-0.678141025604516*err - 0.3205521521864134*lmM - 0.0800290327210927*lmM*lmM - 
           0.010030864197530864*lmM*lmM*lmM) + 
           (-1 + nf)*(-745.7207145811878*err - 358.29765085086774*lmM - 87.39262571554698*lmM*lmM - 
           11.883873456790122*lmM*lmM*lmM - 1.2705439814814814*lmM*lmM*lmM*lmM) + 
           (-1. + nf)*(-1. + nf)*(43.396250117985666*err + 20.528466368867228*lmM + 4.971905254812516*lmM*lmM + 
           0.6304012345679011*lmM*lmM*lmM + 0.06655092592592593*lmM*lmM*lmM*lmM);
     return erg;
}

// Compute \Delta using eq.(14) of [RunDec]
double CRunDec::fDelta(double mOS,double mq[]){
     double erg=0.0;
     for (int i=0; i<4; i++){
       double x= (double) mq[i]/mOS;
       if(x>1.||x<0.){
         cout<<"\\Delta "<<x<<" IS CALLED; THE FUNCTION IS NOT IMPELEMENTED "<<
               "FOR ARGUMENTS OUTSIDE THE INTERVAL [0,1]."<<endl;
         RETURN
       }
       erg+= Pi*Pi*x/8. - 0.597*x*x + 0.230*x*x*x;
     }
     return erg;
}       
// z_m^inv
double CRunDec::fZmInvM(double nl){
     double erg;
     erg=(8481925./93312. + 
       (137.*nl)/216. + (652841.*Pi*Pi)/38880. - (nl*Pi*Pi)/27. - 
       (695.*Pi*Pi*Pi*Pi)/7776. - (575.*Pi*Pi*log(2))/162. - 
       (22.*Pi*Pi*log(2)*log(2))/81. - 
       (55.*log(2)*log(2)*log(2)*log(2))/162. - (220.*A4)/27. - 
       nl*nl*(-2353./23328. - (13.*Pi*Pi)/324. - (7.*Zeta3)/54.) + 
       (58.*Zeta3)/27. - 
       (1439.*Pi*Pi*Zeta3)/432. - nl*(246643./23328. + (967.*Pi*Pi)/648. - 
       (61.*Pi*Pi*Pi*Pi)/1944. + (11.*Pi*Pi*log(2))/81. - 
       (2.*Pi*Pi*log(2)*log(2))/81. - 
       log(2)*log(2)*log(2)*log(2)/81. - (8.*A4)/27. + 
       (241.*Zeta3)/72.) + 
       (1975.*Zeta5)/216.);
     return erg;
}

// Function: double CRunDec::mMS2mOS(double MS, double mq[], double asmu,
//                           double mu,int nl)
double CRunDec::mMS2mOS(double MS, double mq[], double asmu,double mu,int nl, double err){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl == 4 && (Nf<4||Nf>6)){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR 4 LOOPS and "<< nl << " FLAVORS"<<endl;
       RETURN  
     }
     double sum[5];
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fOsFromMs1(mu, MS))/Pi;
     sum[2]=asmu*asmu*((this-> fOsFromMs2(mu, MS, Nf-1))
           +4.*(this->fDelta(MS,mq))/3.)/(Pi*Pi);
     sum[3]=asmu*asmu*asmu*(this-> fOsFromMs3(mu, MS,Nf-1)+
           this->fZmInvM(Nf-1))/(Pi*Pi*Pi);   
     sum[4]=asmu*asmu*asmu*asmu*(this->fOsFromMs4(mu,MS,Nf, err))/(Pi*Pi*Pi*Pi);  
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{ 
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     return MS*erg;
}

// Function: double CRunDec::mMS2mOSmod(double MS, double mq[], double asmu,
//                           double mu,int nl), need for several other mass relations
double CRunDec::mMS2mOSmod(double MS, double mq[], double asmu,double mu,int nf, int nl, double err){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl == 4 && (nf<4||nf>6)){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR 4 LOOPS and "<< nl << " FLAVORS"<<endl;
       RETURN  
     }
     SetConstants(nf);
     double sum[5];
     double lmm = log((mu*mu)/(MS*MS));
     sum[0]=(double) 1.;


     sum[1]=asmu*(this ->fOsFromMs1(mu, MS))/(Pi);
     sum[2]=asmu*asmu*((this-> fOsFromMs2(mu, MS, Nf-1))+4.*(this->fDelta(MS,mq))/3.)/(Pi*Pi) +
            asmu*asmu*(lmm)/(6.*Pi)*(this ->fOsFromMs1(mu, MS))/(Pi);
     sum[3]=asmu*asmu*asmu*(this-> fOsFromMs3(mu, MS,Nf-1)+this->fZmInvM(Nf-1))/(Pi*Pi*Pi)+
            asmu*asmu*asmu*2.*(lmm)/(6.*Pi)*((this-> fOsFromMs2(mu, MS, Nf-1))+4.*(this->fDelta(MS,mq))/3.)/(Pi*Pi)+
            asmu*asmu*asmu*(-11./72. + (11.*lmm)/24. + lmm*lmm/36.)/(Pi*Pi)*(this ->fOsFromMs1(mu, MS))/(Pi);
     sum[4]=asmu*asmu*asmu*asmu*(this->fOsFromMs4(mu,MS,Nf, err))/(Pi*Pi*Pi*Pi) + 
			asmu*asmu*asmu*asmu*3.*(lmm)/(6.*Pi)*(this-> fOsFromMs3(mu, MS,Nf-1)+this->fZmInvM(Nf-1))/(Pi*Pi*Pi) + 
            asmu*asmu*asmu*asmu*2.*(-11./72. + (11.*lmm)/24. + lmm*lmm/36.)/(Pi*Pi)*
                                  ((this-> fOsFromMs2(mu, MS, Nf-1))+4.*(this->fDelta(MS,mq))/3.)/(Pi*Pi) +
            asmu*asmu*asmu*asmu*(-564731./124416. + (2645.*lmm)/1728. + (167.*lmm*lmm)/576. + lmm*lmm*lmm/216. + 
   					(2633./31104. - (67.*lmm)/576. + lmm*lmm/36.)*nl + (82043.*Zeta3)/27648.)/(Pi*Pi*Pi)*(this ->fOsFromMs1(mu, MS))/(Pi) + 
            asmu*asmu*asmu*asmu*(lmm)/(6.*Pi)*(lmm)/(6.*Pi)*((this-> fOsFromMs2(mu, MS, Nf-1))+4.*(this->fDelta(MS,mq))/3.)/(Pi*Pi);  
  
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{ 
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     return MS*erg;
}      

// Coefficients of eq.(16) of [RunDec]
double CRunDec::fMumFromOs1(void){
     return (double) -cf;
}
  
double CRunDec::fMumFromOs2(void){
     double erg;
     erg= (double) ((-1111.*ca*cf)/384. + (199.*cf*cf)/128. + (143.*cf*tr)/96. +
     (71.*cf*(Nf-1)*tr)/96. + (ca*cf*Zeta2)/2. - (15.*cf*cf*Zeta2)/8. - 
     (3.*ca*cf*log(2)*Zeta2)/2. +
     3.*cf*cf*log(2)*Zeta2 - cf*tr*Zeta2 + (cf*(Nf-1)*tr*Zeta2)/2. + 
     (3.*ca*cf*Zeta3)/8. -
     (3.*cf*cf*Zeta3)/4.);     
     return erg;      
}

// z_m^SI
double CRunDec::fMumFromOs3(void){
     double erg;
     erg= (double) -7172965./93312. - 
	 (293.*(Nf-1))/216. - (618281.*Pi*Pi)/38880. - ((Nf-1)*Pi*Pi)/9. + 
     (695.*Pi*Pi*Pi*Pi)/7776. + (623.*Pi*Pi*log(2))/162. + 
     (22.*Pi*Pi*log(2)*log(2))/81. + 
     (55.*log(2)*log(2)*log(2)*log(2))/162. + (220.*A4)/27. + 
     (Nf-1)*(Nf-1)*(-2353./23328. - (13.*Pi*Pi)/324. - (7.*Zeta3)/54.) - 
     (70.*Zeta3)/27. + 
     (1439.*Pi*Pi*Zeta3)/432. + (Nf-1)*(246643./23328. + (967.*Pi*Pi)/648. - 
     (61.*Pi*Pi*Pi*Pi)/1944. + (11.*Pi*Pi*log(2))/81. - 
     (2*Pi*Pi*log(2)*log(2))/81. - 
     log(2)*log(2)*log(2)*log(2)/81. - (8.*A4)/27. + (241.*Zeta3)/72.) - 
     (1975.*Zeta5)/216.;
     return erg;   
}

double CRunDec::fMumFromOs4(double err){
     double erg;
     erg = err*(-3214.227044839041 + 692.4809215366435*(-1. + Nf) - 41.95978562498058*(-1. + Nf)*(-1. + Nf)
                + 0.678141025604516*(-1. + Nf)*(-1. + Nf)*(-1. + Nf));
     return erg;
}

// Function: double CRunDec::mOS2mSI(double mOS, double mq[], double asM, 
//                           int nl)
double CRunDec::mOS2mSI(double mOS, double mq[], double asM, int nl, double err){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl == 4 && (Nf<4||Nf>6)){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR 4 LOOPS and "<< nl << " FLAVORS"<<endl;
       RETURN  
     }

     double sum[5];
     sum[0]=(double) 1.;
     sum[1]=asM*(this ->fMumFromOs1())/Pi;
     sum[2]=asM*asM*(this-> fMumFromOs2()-4.*fDelta(mOS,mq)/3.)/(Pi*Pi);
     sum[3]=asM*asM*asM*(this-> fMumFromOs3())/(Pi*Pi*Pi);     
     sum[4]=asM*asM*asM*asM*(this->fMumFromOs4(1.0)*err)/(Pi*Pi*Pi*Pi);
     

     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
 
     return mOS*erg;       
       
}

// Function: double CRunDec::mOS2mMSrun(double mOS, double mq[], double asmu, 
//                           double mu, int nl)
double CRunDec::mOS2mMSrun(double mOS, double mq[], double asmu, double mu,
			   int nl){
     double asM=0.0;
     asM= this-> AlphasExact(asmu, mu, mOS, nl);
     double mum= this-> mOS2mSI(mOS, mq, asM, nl);
     double asmum= this-> AlphasExact(asmu, mu, mum, nl);
     double newM= this->mMS2mMS(mum, asmum, asmu, nl);
     return newM;       
}

// Function: double CRunDec::mMS2mOSrun(double mMS, double mq[], double asmu, 
//                           double mu, int nl)
double CRunDec::mMS2mOSrun(double mMS, double mq[], double asmu, double mu,
			   int nl){
     double mNeu = mMS2mSI(mMS, asmu, mu, nl);
     double asmNeu = AlphasExact(asmu, mu, mNeu, nl);
     return mMS2mOS(mNeu, mq, asmNeu, mNeu, nl);   
}

// Coefficients of eq.(19) of [RunDec]
double CRunDec::fRiFromMs(double alpha, double nl){
     double sum[4];
     sum[0]= 1.;
     sum[1]= (4.*alpha)/3.;
     sum[2]= alpha*alpha*((1123./72. - (89.*Nf)/144. - (19.*Zeta3)/6.));
     sum[3]= alpha*alpha*alpha*(6663911./41472. - (118325.*Nf)/7776. +
            (4459.*Nf*Nf)/23328. + 
            (4.*(1123./72. - (89.*Nf)/144. - (19.*Zeta3)/6.))/3. - 
            (408007.*Zeta3)/6912. + 
            (617.*Nf*Zeta3)/216. + (Nf*Nf*Zeta3)/54. - (4.*(-995./72. + 
            (89.*Nf)/144. + (19.*Zeta3)/6.))/3. - 
            (5.*Nf*Zeta4)/12. + (185.*Zeta5)/36.);
     double erg=0.0;
     if(nl==0){
       erg=1.;
     }
     else{
       erg=0.;
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     return (erg);        
}

// Function: double CRunDec::mMS2mRI(double mMS, double asmu, int nl)
double CRunDec::mMS2mRI(double mMS, double asmu, int nl){
     if(nl<0||nl>3){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     return (double) mMS*(this->fRiFromMs((asmu/Pi), nl)); 
}

// Coefficients needed for the transformation of mOS to mMSit
double CRunDec::fHelpmOS2mMSit(double MS, double mOS, double mq[], double asmu,
                                double mu,int nl){
     double sum[4];
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fOsFromMs1(mu, MS))/Pi;
     sum[2]=asmu*asmu*((this-> fOsFromMs2(mu, MS, Nf-1))
           +4.*(this->fDelta(mOS,mq))/3.)/(Pi*Pi);
     sum[3]=asmu*asmu*asmu*(this-> fOsFromMs3(mu, MS,Nf-1)+
           +this->fZmInvM(Nf-1))/(Pi*Pi*Pi);
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{ 
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     return erg;
}

// Function: double CRunDec::mOS2mMSit(double mOS, double mq[], double asmu, 
//                          double mu, int nl)
double CRunDec::mOS2mMSit(double mOS, double mq[], double asmu, double mu,
                           int nl){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN
     }
     double epsilonX= 1e-8;
     double x0=mOS-0.1*mOS; 
     double x1=mOS+0.1*mOS;
     double f0= (x0*(this->fHelpmOS2mMSit(x0,mOS, mq, asmu, mu, nl))-mOS);
     double f1= (x1*(this->fHelpmOS2mMSit(x1,mOS, mq, asmu, mu, nl))-mOS);
     if(f0*f1>0){
       cout<<"WARNING: No root can be calculatet!"<<endl;
       RETURN
     }
     double xTest;
     double fTest;
     do{
       xTest= (x0+x1)/2.;
       fTest= (xTest*(this->fHelpmOS2mMSit(xTest,mOS, mq, asmu, mu, nl))-mOS);
       if(f0*fTest<=0){x1= xTest;}
       else {x0= xTest;}
     }
     while(abs(x1-x0)>= epsilonX);
       double mNeu=xTest;
     return mNeu;               
}

double CRunDec::PSdelta(double asmu, double muf, double mu, int nl, int nloops) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     double lf = log(mu*mu/(muf*muf));
     double api = asmu;
     double ret[5];
     ret[0] = 0.0;
     ret[1] = api*4./3.;
     ret[2] = api*api*(97./9. + nl*(-22./27. - (2.*lf)/9.) + (11.*lf)/3.);
     ret[3] = api*api*api*(33623./216. + 3.*Pi*Pi - (3.*Pi*Pi*Pi*Pi)/16. + (610.*lf)/9.
              + (121.*lf*lf)/12. + nl*nl*(157./243. + (22.*lf)/81. + lf*lf/27.) +
              nl*(-7145./324. - (493.*lf)/54. - (11.*lf*lf)/9. - (13.*Zeta3)/9.) + (11.*Zeta3)/2.);
     ret[4] = api*api*api*api*(3567.723056629293 + (7271.*lf*lf)/24. + (1331.*lf*lf*lf)/48.
              + nl*nl*nl*(-2951./4374. - (157.*lf)/486. - (11.*lf*lf)/162. - 
              lf*lf*lf/162.) + nl*(-701.2303148875468 - (8485.*lf*lf)/144. - (121.*lf*lf*lf)/24.
              + lf*(-253189./864. - (3.*Pi*Pi)/2. + (3.*Pi*Pi*Pi*Pi)/32.
              - (44.*Zeta3)/3.)) + nl*nl*(1751971./46656. + Pi*Pi*Pi*Pi/135. + (773.*lf*lf)/216. +
              (11.*lf*lf*lf)/36. + lf*(15355./864. + (13.*Zeta3)/18.) + 
              (259.*Zeta3)/108.) + lf*(26125./18. + (99.*Pi*Pi)/4. - (99.*Pi*Pi*Pi*Pi)/64. + (363.*Zeta3)/8.));

     double delta = 0.0;
     for(int i = 0; i <= nloops; i++) {
       delta += ret[i];
     }
	 return delta;
}

// Function: double CRunDec::mOS2mPS(double mOS, double mq[], double asmu, double mu, double muf, int nl, int nloops)
double CRunDec::mOS2mPS(double mOS, double mq[], double asmu, double mu, double muf, int nl, int nloops) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }

	 return mOS - muf*PSdelta(asmu/Pi, muf, mu, nl, nloops);
}

// Function: double CRunDec::mMS2mPS(double mMS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err)
double CRunDec::mMS2mPS(double mMS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }

     double delmuf = PSdelta(asmu/Pi, muf, mu, nl, nloops); 
     double exmOS = mMS2mOSmod(mMS, this->mq, asmu, mu, nl+1  , nloops, err);
     return (exmOS - muf*delmuf);
}

// Function: double CRunDec::mPS2mMS(double mPS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err)
double CRunDec::mPS2mMS(double mPS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }

     double lowbound = mPS - mPS/5.;
     double highbound = mPS + mPS/5.;
     bool found = false;
     double f1 = mMS2mPS(lowbound, mq, asmu, mu, muf, nl, nloops, err) - mPS;
     double f2 = mMS2mPS(highbound, mq, asmu, mu, muf, nl, nloops, err) - mPS;
     for(int j = 0; j < 1000; j++) {
       if(f1*f2 < 0.0) {
         found = true;
         break;
       }
       if(fabs(f1) < fabs(f2))
         f1 = mMS2mPS(lowbound += 1.5*(lowbound - highbound), mq, asmu, mu, muf, nl, nloops, err) - mPS;
       else
         f2 = mMS2mPS(lowbound += 1.5*(lowbound - highbound), mq, asmu, mu, muf, nl, nloops, err) - mPS;
     }
     if(found) {
     double acc = 10e-10;
     double dx;
     double xmid;
     double mMS = f1 < 0.0 ? (dx=highbound-lowbound,lowbound) : (dx=lowbound-highbound,highbound);
     for(int j = 0; j < 1000; j++) {
       f2 = mMS2mPS(xmid = mMS+(dx *= 0.5), mq, asmu, mu, muf, nl, nloops, err) - mPS;
       if(f2 <= 0.0)
         mMS = xmid;
       if(fabs(dx) < acc || f2 == 0.0)
         return mMS;
     }
     }
     return 0.0;
}

// Function: double CRunDec::mPS2mSI(double mPS, double mq[], double (*as)(double), double muf, int nl, int nloops, double err)
// The function pointer passed should contain the adress of a function computing alpha_s in dependence of mu
double CRunDec::mPS2mSI(double mPS, double mq[], double (*as)(double), double muf, int nl, int nloops, double err) {
    if(as == NULL) {
      cout << "Pointer to as == NULL! Aborting..." << endl;
      RETURN
    }
	double mMS1 = 0;
	double mMS = mPS;
	double acc = 10e-6; 
	while(fabs(mMS1 - mMS) > acc) {
		mMS1 = mMS; 
        mMS = mPS2mMS(mPS, mq, as(mMS1), mMS1, muf, nl, nloops, err);
	}
	return mMS;
}


double CRunDec::E1p(double mOS, double asmu, double mu, int nl, int nloops) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     double lmm = log(mu/mOS);
     double lmm34 = log(3.*mu/(4.*asmu*mOS));
     double api = asmu/(Pi);
     double ret[5];

     ret[0] = 0.0;
     ret[1] = 1.0;
     ret[2] = api*(97./6. + nl*(-11./9. - (2.*lmm34)/
                 3.) + 11.*lmm34);
     ret[3] = api*api*(1793./12. + (2917./216. - (11.*nl)/18. + nl*nl/54.)*Pi*Pi - 
              (9.*Pi*Pi*Pi*Pi)/32. + (927.*lmm34)/4. + 
              (363.*lmm34*lmm34)/4. + nl*(-1693./72. - 
                (193.*lmm34)/6. - 11.*lmm34*lmm34 - (19.*Zeta3)/2.) + nl*nl*(77./108. + 
                lmm34 + lmm34*lmm34/3. + 
                (2.*Zeta3)/9.) + (275.*Zeta3)/4.);
     ret[4] = api*api*api*(1267919./1728. + Pi*Pi*Pi*Pi*(-723119./51840. + (11.*nl*nl)/
                 1080. - nl*nl*nl/4860. + nl*(59677./77760. + (3.*lmm34)/8.) - (99.*lmm34)/16.) + 
              (4521.*lmm34*lmm34)/2. + (1331.*
                lmm34*lmm34*lmm34)/2. + (114917.*Zeta3)/48. + 
              lmm34*(247675./96. + (3025.*Zeta3)/2.) + 
              Pi*Pi*(265.389067842508 + (865.*lmm)/18. + 
                (26897.*lmm34)/108. + nl*nl*(905./432. + 
                  (11.*lmm34)/9. - (11.*Zeta3)/6.) + 
                nl*nl*nl*(-19./486. - (2.*lmm34)/81. + 
                  Zeta3/27.) + nl*(-397591./7776. - (5095.*lmm34)/162. + (121.*Zeta3)/4.)) + (13432.614375 - 
                3289.906669391583*nl - (1000.*nl*nl*nl)/729. + nl*nl*
                 ((14002./81. - (416.*Zeta3)/3.)/3. + (3.*(12541./243. + (64.*Pi*Pi*Pi*Pi)/
                      135. + (368.*Zeta3)/3.))/4.))/32. + nl*(-52033./288. - 
                (10955.*lmm34*lmm34)/24. - 
                121.*lmm34*lmm34*lmm34 + lmm34*
                 (-166309./288. - (902.*Zeta3)/3.) - (8797.*Zeta3)/18. - 
                363.*Zeta5) + nl*nl*nl*(-98./729. - (5.*lmm34*lmm34)/
                 9. - (4.*lmm34*lmm34*lmm34)/27. + 
                lmm34*(-50./81. - (8.*Zeta3)/27.) - 
                (44.*Zeta3)/81. - (4.*Zeta5)/9.) + (3993.*Zeta5)/2. + 
              nl*nl*(3073./288. + (1027.*lmm34*lmm34)/36. + 
                (22.*lmm34*lmm34*lmm34)/3. + (3239.*Zeta3)/108. + 
                lmm34*(10351./288. + (158.*Zeta3)/9.) + 
                22.*Zeta5));


     double E = 0.0;
     for(int i = 0; i <= nloops; i++) {
       E += ret[i];
     }
	 return -E*(4.*asmu*asmu*mOS)/9.;
}

// Function: double CRunDec::mOS2m1S(double mOS, double mq[], double asmu, double mu, int nl, int nloops)
double CRunDec::mOS2m1S(double mOS, double mq[], double asmu, double mu, int nl, int nloops) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
	 return mOS + 0.5*E1p(mOS, asmu, mu, nl, nloops);
}

// Function: double CRunDec::mMS2m1S(double mMS, double mq[], double asmu, double mu, int nl, int nloops, double err)
// Derived from RunDec
double CRunDec::mMS2m1S(double mMS, double mq[], double asmu, double mu, int nl, int nloops, double err) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     if(nl == 4 && (nl<4||nl>6)){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR 4 LOOPS and "<< nl << " FLAVORS"<<endl;
       RETURN  
     }
     double lmm = log(mu*mu/(mMS*mMS));
     double log34 = log((3.*mu)/(4.*asmu*mMS));
     double m1S = mMS;
     switch(nloops) {
       case 1: m1S += asmu*mMS*(12. - 2.*asmu*Pi + 9.*lmm)/(9.*Pi);
               break;
       case 2: m1S += asmu*mMS*(-(2.*asmu)/9. + asmu*asmu*(-291. + 22.*nl - 198.*log34 + 
                      12.*nl*log34 - 3.*(8. + 6.*lmm))/(81.*Pi) + (4. + 3.*lmm)/(3.*Pi) +
                      (asmu*(2763. - 142.*nl + 96.*Pi*Pi - 16.*nl*Pi*Pi + 32.*Pi*Pi*ln2 + 
                      2036.*lmm - 104.*nl*lmm + 564.*lmm*lmm - 24.*nl*lmm*lmm - 48.*Zeta3))/(288.*Pi*Pi));
               break;
       case 3: m1S += asmu*mMS*(-(2.*asmu)/9. + asmu*asmu*(-291. + 22.*nl -  198.*log34 + 
                      12.*nl*log34)/(81.*Pi) + (4. + 3.*lmm)/(3.*Pi) - (2.*asmu*asmu*(4. + 3.*lmm))/(27.*Pi) + 
                      asmu*asmu*asmu*(-372. + 40.*nl - 792.*log34 + 
                      48.*nl*log34 - 279.*lmm + 30.*nl*lmm - 
                      594.*log34*lmm + 36.*nl*log34*lmm)/(243.*Pi*Pi)
                      + (asmu*(2763. - 142.*nl + 96.*Pi*Pi - 16.*nl*Pi*Pi + 32.*Pi*Pi*ln2 + 2036.*lmm - 104.*nl*lmm + 
                      564.*lmm*lmm - 24.*nl*lmm*lmm - 48.*Zeta3))/(288.*Pi*Pi) + asmu*asmu*asmu*(-2763.+ 142.*nl - 
                      96.*Pi*Pi + 16.*nl*Pi*Pi - 32.*Pi*Pi*ln2 -
                      2036.*lmm + 104.*nl*lmm - 564.*lmm*lmm + 
                      24.*nl*lmm*lmm + 48.*Zeta3)/(1296.*Pi*Pi)
                      + asmu*asmu*asmu*(-129096. + 20316.*nl - 616.*nl*nl - 
                      11668.*Pi*Pi + 528.*nl*Pi*Pi - 16.*nl*nl*Pi*Pi + 243.*Pi*Pi*Pi*Pi - 
                      200232.*log34 + 27792.*nl*log34 - 864.*nl*nl*log34 - 
                      78408.*log34*log34 + 9504.*nl*log34*log34 - 
                      288.*nl*nl*log34*log34 - 59400.*Zeta3 + 8208.*nl*Zeta3 - 
                      192.*nl*nl*Zeta3)/(3888*Pi*Pi) + 
                      (asmu*asmu*(42314585. - 4636940.*nl + 47060.*nl*nl + 7834092.*Pi*Pi - 
                      713520.*nl*Pi*Pi + 18720.*nl*nl*Pi*Pi - 41700.*Pi*Pi*Pi*Pi + 14640.*nl*Pi*Pi*Pi*Pi - 
                      1656000.*Pi*Pi*ln2 - 63360.*nl*Pi*Pi*ln2 - 126720.*Pi*Pi*ln2*ln2 + 
                      11520.*nl*Pi*Pi*ln2*ln2 - 158400.*ln2*ln2*ln2*ln2 + 5760.*nl*ln2*ln2*ln2*ln2 + 
                      33620760.*lmm - 3723120.*nl*lmm + 64080.*nl*nl*lmm + 1010880.*Pi*Pi*lmm - 220320.*nl*Pi*Pi*lmm + 8640.*nl*nl*Pi*Pi*lmm + 
                      336960.*Pi*Pi*ln2*lmm - 17280.*nl*Pi*Pi*ln2*lmm + 11726100.*lmm*lmm - 1247400.*nl*lmm*lmm + 28080.*nl*nl*lmm*lmm + 
                      2009880.*lmm*lmm*lmm - 185760.*nl*lmm*lmm*lmm + 4320.*nl*nl*lmm*lmm*lmm - 3801600.*A4 + 138240.*nl*A4 + 1002240.*Zeta3 - 
                      1561680.*nl*Zeta3 + 60480.*nl*nl*Zeta3 - 1554120.*Pi*Pi*Zeta3 - 894240.*lmm*Zeta3 - 362880.*nl*lmm*Zeta3 + 4266000.*Zeta5))/
                      (466560.*Pi*Pi*Pi));
                      break;
       case 4: m1S += - (2*asmu*asmu*mMS)/9 + (-291*asmu*asmu*asmu*mMS + 22*asmu*asmu*asmu*mMS*nl - 
                      198*asmu*asmu*asmu*mMS*log34 + 
                      12*asmu*asmu*asmu*mMS*nl*log34)/(81*Pi) + 
                      (asmu*mMS*(4 + 3*lmm))/(3*Pi) - 
                      (2*(4*asmu*asmu*asmu*mMS + 3*asmu*asmu*asmu*mMS*lmm))/(27*Pi) + 
                      (-372*asmu*asmu*asmu*asmu*mMS + 40*asmu*asmu*asmu*asmu*mMS*nl - 
                      792*asmu*asmu*asmu*asmu*mMS*log34 + 
                      48*asmu*asmu*asmu*asmu*mMS*nl*log34 - 
                      279*asmu*asmu*asmu*asmu*mMS*lmm + 30*asmu*asmu*asmu*asmu*mMS*nl*lmm - 
                      594*asmu*asmu*asmu*asmu*mMS*log34*lmm + 
                      36*asmu*asmu*asmu*asmu*mMS*nl*log34*lmm)/(243*Pi*Pi)
                      + 1.5139171475883256E-8*(-3.6643425048120157E6*asmu*asmu*asmu*asmu*mMS - 
                      3.670703784417217E9*asmu*asmu*asmu*asmu*asmu*mMS + 292307.7010004199*asmu*asmu*asmu*asmu*mMS*nl + 
                      6.48095550442438E8*asmu*asmu*asmu*asmu*asmu*mMS*nl - 
                      3.6315792483681194E7*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl + 
                      595364.9557740748*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*nl - 
                      2.245324461595123E8*asmu*asmu*asmu*asmu*asmu*mMS*log(mu/mMS) - 
                      3.062101652047988E9*asmu*asmu*asmu*asmu*asmu*mMS*log34 + 
                      5.902121099662828E8*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34 - 
                      3.325469404123189E7*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34 + 
                      576205.4055744545*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*nl*log34 - 
                      1.1274215864457371E9*asmu*asmu*asmu*asmu*asmu*mMS*log34*log34 + 
                      2.230344299450968E8*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34*log34 - 
                      1.3715687063042711E7*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*log34 + 
                      263004.5457918066*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*nl*log34*log34 - 
                      3.1505314540400505E8*asmu*asmu*asmu*asmu*asmu*mMS*log34*log34*log34 + 
                      5.728239007345548E7*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34*log34*log34 - 
                      3.4716600044518467E6*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*log34*log34 + 
                      70134.54554448176*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*nl*log34*log34*log34 + 
                      7.228435160639514E7*asmu*asmu*asmu*asmu*mMS*lmm - 
                      1.061086986213349E8*asmu*asmu*asmu*asmu*asmu*mMS*lmm - 
                      9.657775840868242E6*asmu*asmu*asmu*asmu*mMS*nl*lmm + 
                      1.3033590170779591E7*asmu*asmu*asmu*asmu*asmu*mMS*nl*lmm + 
                      221297.1713327142*asmu*asmu*asmu*asmu*mMS*nl*nl*lmm - 
                      323617.3768070282*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*lmm - 
                      6.0602822464077026E7*asmu*asmu*asmu*asmu*asmu*mMS*log34*
                      lmm + 8.924620920535302E6*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34*lmm - 
                      271771.36398486677*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*
                      lmm - 4.2961792555091605E7*asmu*asmu*asmu*asmu*asmu*mMS*log34*log34*lmm + 
                      5.20749000667777E6*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34*log34*
                      lmm - 157802.72747508393*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*log34*lmm + 
                      3.6846948160874695E7*asmu*asmu*asmu*asmu*mMS*lmm*lmm - 
                      1.4084441353288308E7*asmu*asmu*asmu*asmu*asmu*mMS*lmm*lmm - 
                      4.890850849173104E6*asmu*asmu*asmu*asmu*mMS*nl*lmm*lmm + 
                      1.826785740978923E6*asmu*asmu*asmu*asmu*asmu*mMS*nl*lmm*lmm + 
                      108536.53795824476*asmu*asmu*asmu*asmu*mMS*nl*nl*lmm*lmm - 
                      50409.20461009625*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*lmm*lmm - 
                      1.0198001263077298E7*asmu*asmu*asmu*asmu*asmu*mMS*log34*
                      lmm*lmm + 1.0520181831672261E6*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34*lmm*lmm - 
                      26300.454579180656*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*
                      lmm*lmm + 9.281070040868293E6*asmu*asmu*asmu*asmu*mMS*
                      lmm*lmm*lmm - 2.039381082160633E6*asmu*asmu*asmu*asmu*asmu*mMS*
                      lmm*lmm*lmm - 919055.1091884745*asmu*asmu*asmu*asmu*mMS*nl*
                      lmm*lmm*lmm + 188486.5911507947*asmu*asmu*asmu*asmu*asmu*mMS*nl*
                      lmm*lmm*lmm + 20406.00584022232*asmu*asmu*asmu*asmu*mMS*nl*nl*
                      lmm*lmm*lmm - 4383.409096530109*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*
                      lmm*lmm*lmm + 1.3530751564824337E6*asmu*asmu*asmu*asmu*mMS*
                      lmm*lmm*lmm*lmm - 130284.49882603479*asmu*asmu*asmu*asmu*mMS*nl*
                      lmm*lmm*lmm*lmm + 3139.3855138803565*asmu*asmu*asmu*asmu*mMS*nl*nl*
                      lmm*lmm*lmm*lmm + 678107.270998157*asmu*asmu*asmu*asmu*mMS*
                      this->fOsFromMs4(mu,mMS,nl+1,err)) + 
                      (asmu*asmu*mMS*(2763 - 142*nl + 96*Pi*Pi - 16*nl*Pi*Pi + 32*Pi*Pi*ln2 + 
                      2036*lmm - 104*nl*lmm + 564*lmm*lmm - 24*nl*lmm*lmm - 48*Zeta3))/
                      (288*Pi*Pi) + (-2763*asmu*asmu*asmu*asmu*mMS + 142*asmu*asmu*asmu*asmu*mMS*nl - 
                      96*asmu*asmu*asmu*asmu*mMS*Pi*Pi + 16*asmu*asmu*asmu*asmu*mMS*nl*Pi*Pi - 
                      32*asmu*asmu*asmu*asmu*mMS*Pi*Pi*ln2 - 2036*asmu*asmu*asmu*asmu*mMS*lmm + 
                      104*asmu*asmu*asmu*asmu*mMS*nl*lmm - 564*asmu*asmu*asmu*asmu*mMS*lmm*lmm + 
                      24*asmu*asmu*asmu*asmu*mMS*nl*lmm*lmm + 48*asmu*asmu*asmu*asmu*mMS*Zeta3)/(1296*Pi*Pi)
                      + (-129096*asmu*asmu*asmu*asmu*mMS + 20316*asmu*asmu*asmu*asmu*mMS*nl - 616*asmu*asmu*asmu*asmu*mMS*nl*nl - 
                      11668*asmu*asmu*asmu*asmu*mMS*Pi*Pi + 528*asmu*asmu*asmu*asmu*mMS*nl*Pi*Pi - 
                      16*asmu*asmu*asmu*asmu*mMS*nl*nl*Pi*Pi + 243*asmu*asmu*asmu*asmu*mMS*Pi*Pi*Pi*Pi - 
                      200232*asmu*asmu*asmu*asmu*mMS*log34 + 
                      27792*asmu*asmu*asmu*asmu*mMS*nl*log34 - 
                      864*asmu*asmu*asmu*asmu*mMS*nl*nl*log34 - 
                      78408*asmu*asmu*asmu*asmu*mMS*log34*log34 + 
                      9504*asmu*asmu*asmu*asmu*mMS*nl*log34*log34 - 
                      288*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*log34 - 
                      59400*asmu*asmu*asmu*asmu*mMS*Zeta3 + 8208*asmu*asmu*asmu*asmu*mMS*nl*Zeta3 - 
                      192*asmu*asmu*asmu*asmu*mMS*nl*nl*Zeta3)/(3888*Pi*Pi) + 
                      (asmu*asmu*asmu*mMS*(42314585 - 4636940*nl + 47060*nl*nl + 7834092*Pi*Pi - 
                      713520*nl*Pi*Pi + 18720*nl*nl*Pi*Pi - 41700*Pi*Pi*Pi*Pi + 14640*nl*Pi*Pi*Pi*Pi - 
                      1656000*Pi*Pi*ln2 - 63360*nl*Pi*Pi*ln2 - 126720*Pi*Pi*ln2*ln2 + 
                      11520*nl*Pi*Pi*ln2*ln2 - 158400*ln2*ln2*ln2*ln2 + 5760*nl*ln2*ln2*ln2*ln2 + 
                      33620760*lmm - 3723120*nl*lmm + 
                      64080*nl*nl*lmm + 1010880*Pi*Pi*lmm - 
                      220320*nl*Pi*Pi*lmm + 8640*nl*nl*Pi*Pi*lmm + 
                      336960*Pi*Pi*ln2*lmm - 
                      17280*nl*Pi*Pi*ln2*lmm + 11726100*lmm*lmm - 
                      1247400*nl*lmm*lmm + 28080*nl*nl*lmm*lmm + 
                      2009880*lmm*lmm*lmm - 185760*nl*lmm*lmm*lmm + 
                      4320*nl*nl*lmm*lmm*lmm - 3801600*A4 + 
                      138240*nl*A4 + 1002240*Zeta3 - 1561680*nl*Zeta3 + 
                      60480*nl*nl*Zeta3 - 1554120*Pi*Pi*Zeta3 - 
                      894240*lmm*Zeta3 - 362880*nl*lmm*Zeta3 + 
                      4266000*Zeta5))/(466560*Pi*Pi*Pi);
                      break;
     }



     return m1S;  
}

// Function: double CRunDec::m1S2mMS(double m1S, double mq[], double asmu, double mu, int nl, int nloops, double err)
double CRunDec::m1S2mMS(double m1S, double mq[], double asmu, double mu, int nl, int nloops, double err) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }

     double lowbound = m1S - m1S/5.;
     double highbound = m1S + m1S/5.;
     bool found = false;
     double f1 = mMS2m1S(lowbound, mq, asmu, mu, nl, nloops, err) - m1S;
     double f2 = mMS2m1S(highbound,mq, asmu, mu, nl, nloops, err) - m1S;
     for(int j = 0; j < 1000; j++) {
       if(f1*f2 < 0.0) {
         found = true;
         break;
       }
       if(fabs(f1) < fabs(f2))
         f1 = mMS2m1S(lowbound += 1.5*(lowbound - highbound), mq, asmu, mu, nl, nloops, err) - m1S;
       else
         f2 = mMS2m1S(lowbound += 1.5*(lowbound - highbound), mq, asmu, mu, nl, nloops, err) - m1S;
     }
     if(found) {
     double acc = 10e-10;
     double dx;
     double xmid;
     double mMS = f1 < 0.0 ? (dx=highbound-lowbound,lowbound) : (dx=lowbound-highbound,highbound);
     for(int j = 0; j < 1000; j++) {
       f2 = mMS2m1S(xmid = mMS+(dx *= 0.5), mq, asmu, mu, nl, nloops, err) - m1S;
       if(f2 <= 0.0)
         mMS = xmid;
       if(fabs(dx) < acc || f2 == 0.0)
         return mMS;
     }
     }
     return 0.0;
}

// Function: double CRunDec::m1S2mSI(double m1S, double mq[], double (*as)(double), int nl, int nloops, double err)
// The function pointer passed should contain the adress of a function computing alpha_s in dependence of mu
double CRunDec::m1S2mSI(double m1S, double mq[], double (*as)(double), int nl, int nloops, double err) {
    if(as == NULL) {
      cout << "Pointer to as == NULL! Aborting..." << endl;
      RETURN
    }
	double mMS1 = 0;
	double mMS = m1S;
	double acc = 10e-6; 
	while(fabs(mMS1 - mMS) > acc) {
		mMS1 = mMS; 
        mMS = m1S2mMS(m1S, mq, as(mMS1), mMS1, nl, nloops, err);
	}
	return mMS;
}

double CRunDec::exOS2RS(double api, double mmu, double nnuf, int nnl, int nloops) {
	if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     double lmm = log(mmu*mmu/(nnuf*nnuf));
     double ret[5];
     double Nm[] = {1.0, 1.0, 1.0, 0.563, 0.547, 0.527};
     double nu[] = {1.0, 1.0, 1.0, 0.3951, 0.3696, 0.3289};
     double ctil[4][6];
     ctil[1][3] = -0.1638;
     ctil[1][4] = -0.1054;
     ctil[1][5] = 0.0238;
     ctil[2][3] = 0.2372;
     ctil[2][4] = 0.2736;
     ctil[2][5] = 0.3265;
     ctil[3][3] = -0.1205;
     ctil[3][4] = -0.161;
     ctil[3][5] = -0.2681;
     ctil[0][0] = 1;
     ctil[0][1] = 1;
     ctil[0][2] = 1;
     ctil[0][3] = 1;
     ctil[0][4] = 1;
     ctil[0][5] = 1;


     ret[0] = 0.0;
     ret[1] = api*nnuf*Pi*(1. + ctil[1][nnl] + ctil[2][nnl] + ctil[3][nnl])*Nm[nnl];
     ret[2] = api*api*(((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/
        2. - (nnuf*Pi*(1. + ctil[1][nnl] + ctil[2][nnl] + ctil[3][nnl])*
         (-33.*lmm + 2.*nnl*lmm)*Nm[nnl])/12.);

     ret[3] = api*api*api*(((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(3. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/
        4. + nnuf*Pi*(1. + ctil[1][nnl] + ctil[2][nnl] + ctil[3][nnl])*
        lmm*((102. - (38.*nnl)/3.)/16. + 
         ((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*lmm)/16.)*Nm[nnl] + 
       ((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*
         ((11.*lmm)/2. - (nnl*lmm)/3.)*Nm[nnl])/2.);

     ret[4] = api*api*api*api*(((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(3. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(4. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/
        8. + ((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(3. + nu[nnl])/tgamma(1. + nu[nnl]))*
         ((33.*lmm)/4. - (nnl*lmm)/2.)*Nm[nnl])/4. + 
       nnuf*Pi*(1. + ctil[1][nnl] + ctil[2][nnl] + ctil[3][nnl])*
        lmm*((2857./2. - (5033.*nnl)/18. + (325.*nnl*nnl)/54.)/64. + 
         (5.*(102. - (38.*nnl)/3.)*(11. - (2.*nnl)/3.)*lmm)/128. + 
         ((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*lmm*lmm)/64.)*Nm[nnl] + 
       ((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*
         (612.*lmm - 76.*nnl*lmm + 
          1089.*lmm*lmm - 132.*nnl*lmm*lmm + 
          4.*nnl*nnl*lmm*lmm)*Nm[nnl])/96.);


     double res = 0.0;
     for(int i = 0; i <= nloops; i++) {
       res += ret[i];
     }
	 return res;
}

double CRunDec::exOS2RSp(double api, double mmu, double nnuf, int nnl, int nloops) {
	if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     double lmm = log(mmu*mmu/(nnuf*nnuf));
     double ret[5];
     double Nm[] = {1.0, 1.0, 1.0, 0.563, 0.547, 0.527};
     double nu[] = {1.0, 1.0, 1.0, 0.3951, 0.3696, 0.3289};
     double ctil[4][6];
     ctil[1][3] = -0.1638;
     ctil[1][4] = -0.1054;
     ctil[1][5] = 0.0238;
     ctil[2][3] = 0.2372;
     ctil[2][4] = 0.2736;
     ctil[2][5] = 0.3265;
     ctil[3][3] = -0.1205;
     ctil[3][4] = -0.161;
     ctil[3][5] = -0.2681;
     ctil[0][0] = 1;
     ctil[0][1] = 1;
     ctil[0][2] = 1;
     ctil[0][3] = 1;
     ctil[0][4] = 1;
     ctil[0][5] = 1;


     ret[0] = 0.0;
     ret[1] = 0.0;
     ret[2] = api*api*((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
         tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
         tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
         tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/2.;

     ret[3] = api*api*api*(((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(3. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/
        4. + ((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*
         ((11.*lmm)/2. - (nnl*lmm)/3.)*Nm[nnl])/2.);

     ret[4] = api*api*api*api*(((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(3. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(4. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/
        8. + ((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(3. + nu[nnl])/tgamma(1. + nu[nnl]))*
         ((33.*lmm)/4. - (nnl*lmm)/2.)*Nm[nnl])/4. + 
       ((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*
         (612.*lmm - 76.*nnl*lmm + 
          1089.*lmm*lmm - 132.*nnl*lmm*lmm + 
          4.*nnl*nnl*lmm*lmm)*Nm[nnl])/96.);

     double res = 0.0;
     for(int i = 0; i <= nloops; i++) {
       res += ret[i];
     }
	 return res;
}

// Function: CRunDec::mOS2mRS(double mOS, double mq[], double asmu, double mu, double nuf, int nl, int nloops, bool prime)
double CRunDec::mOS2mRS(double mOS, double mq[], double asmu, double mu, double nuf, int nl, int nloops, bool prime) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     if(prime) {
	   return mOS - exOS2RSp(asmu/Pi, mu, nuf, nl, nloops);
     } else {
	   return mOS - exOS2RS(asmu/Pi, mu, nuf, nl, nloops);
     }
}

// Function: double CRunDec::mMS2mRS(double mMS, double mq[], double asmu, double mu, double nuf, int nl, int nloops, double err, bool prime)
double CRunDec::mMS2mRS(double mMS, double mq[], double asmu, double mu, double nuf, int nl, int nloops, double err, bool prime) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     if(prime) {
	   return mMS2mOSmod(mMS, mq, asmu, mu, nl+1, nloops, err) - exOS2RSp(asmu/Pi, mu, nuf, nl, nloops);
     } else {
	   return mMS2mOSmod(mMS, mq, asmu, mu, nl+1, nloops, err) - exOS2RS(asmu/Pi, mu, nuf, nl, nloops);
     }
}

// Function: double CRunDec::mRS2mMS(double mRS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err, bool prime)
double CRunDec::mRS2mMS(double mRS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err, bool prime) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }

     double lowbound = mRS - mRS/5.;
     double highbound = mRS + mRS/5.;
     bool found = false;
     double f1 = mMS2mRS(lowbound, mq, asmu, mu, muf, nl, nloops, err, prime) - mRS;
     double f2 = mMS2mRS(highbound, mq, asmu, mu, muf, nl, nloops, err, prime) - mRS;
     for(int j = 0; j < 1000; j++) {
       if(f1*f2 < 0.0) {
         found = true;
         break;
       }
       if(fabs(f1) < fabs(f2))
         f1 = mMS2mRS(lowbound += 1.5*(lowbound - highbound), mq, asmu, mu, muf, nl, nloops, err, prime) - mRS;
       else
         f2 = mMS2mRS(lowbound += 1.5*(lowbound - highbound), mq, asmu, mu, muf, nl, nloops, err, prime) - mRS;
     }
     if(found) {
     double acc = 10e-10;
     double dx;
     double xmid;
     double mMS = f1 < 0.0 ? (dx=highbound-lowbound,lowbound) : (dx=lowbound-highbound,highbound);
     for(int j = 0; j < 1000; j++) {
       f2 = mMS2mRS(xmid = mMS+(dx *= 0.5), mq, asmu, mu, muf, nl, nloops, err, prime) - mRS;
       if(f2 <= 0.0)
         mMS = xmid;
       if(fabs(dx) < acc || f2 == 0.0)
         return mMS;
     }
     }
     return 0.0;
}

// Function: double CRunDec::mRS2mSI(double mRS, double mq[], double (*as)(double), double muf, int nl, int nloops, double err, bool prime)
// The function pointer passed should contain the adress of a function computing alpha_s in dependence of mu
double CRunDec::mRS2mSI(double mRS, double mq[], double (*as)(double), double muf, int nl, int nloops, double err, bool prime) {
    if(as == NULL) {
      cout << "Pointer to as == NULL! Aborting..." << endl;
      RETURN
    }
	double mMS1 = 0;
	double mMS = mRS;
	double acc = 10e-6; 
	while(fabs(mMS1 - mMS) > acc) {
		mMS1 = mMS; 
        mMS = mRS2mMS(mRS, mq, as(mMS1), mMS1, muf, nl, nloops, err, prime);
	}
	return mMS;
}

// Function: double CRunDec::mMS2mRGImod(double mMS, double asmu, int nl)
// See 'mMS2mRGI' but 'Alphas/Pi -> AlphaS*2*Beta0/Pi'
double CRunDec::mMS2mRGImod(double mMS, double asmu, int nl){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl==0){
       return (double) mMS;
     }
     double bet0= 11./4. - Nf/6.;     
     double cAsmu= this->fSetcx((2*bet0*asmu)/Pi, nl); 
     return (double) mMS/cAsmu;     
}

// Coefficients of eq.(22) of [RunDec]
double CRunDec::fas5to6os(double A,double mass,double mu,double nlq,double nl){
     double log2 = log(2.);
     double lmM  = log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=A*(-lmM)/6.;
     sum[2]=A*A*(-7./24. - (19.*lmM)/24. + (lmM*lmM)/36.);
     sum[3]=A*A*A*(-58933./124416. - (8521.*lmM)/1728. - (131.*lmM*lmM)/576. -
                     (lmM*lmM*lmM)/216. +
       nlq*(2479./31104. + (409.*lmM)/1728. + Zeta2/9.) - (2.*Zeta2)/3. - 
       (2.*Zeta2*log(2))/9. -
       (80507.*Zeta3)/27648.);
     sum[4]=
       + pow(A,4)*(
		   -39.499192691910871413 + (2965705*A4)/54432. - (121*A5)/36. - 
       (7693*pow(lmM,2))/1152. - (8371*pow(lmM,3))/10368. + 
       pow(lmM,4)/1296. - (697121*Zeta2)/19440. - (49*log2*Zeta2)/54. - 
       (243892631*Zeta4)/8.70912e6 + (605*Zeta2*Zeta4)/2688. + 
       (330575*Zeta5)/41472. + (587*Zeta2*log2)/81. - 
       (2057*Zeta4*log2)/576. - (2699593*Zeta2*pow(log2,2))/217728. - 
       (121*Zeta2*pow(log2,3))/432. + 
       lmM*(-26.38581988383059 - (29*Zeta2)/9. - (29*log2*Zeta2)/27. - 
          (2439119*Zeta3)/165888.) + 
       pow(nlq,2)*(-0.09432401513203018 - (493*pow(lmM,2))/20736. + 
          lmM*(-0.008996699245541839 - Zeta2/27.) - (13*Zeta2)/162. - 
          (19*Zeta3)/1728.) + (1439*Zeta2*Zeta3)/216. + 
       nlq*(2.375194240826475 + (173*A4)/5184. + (6661*pow(lmM,2))/10368. + 
          (107*pow(lmM,3))/1728. + (557*Zeta2)/162. - 
          (697709*Zeta4)/165888. + (115*Zeta5)/576. + (22*Zeta2*log2)/81. - 
	    (1709*Zeta2*pow(log2,2))/20736. + 
          (173*pow(log2,4))/124416. + (4756441*Zeta3)/995328. + 
          lmM*(2.975080911351166 + (41*Zeta2)/54. + (2*log2*Zeta2)/27. + 
             (132283*Zeta3)/82944.))
		   );

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsDownOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsDownOS(double als, double massth, double muth, int nl){
       if(nl<1||nl>5){
         cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
         RETURN 
       }
       double erg=(this->fas5to6os(als/Pi, massth, muth, Nf,nl));
       return als*erg;
}

// Coefficients of eq.(25) of [RunDec]
double CRunDec::fas6to5os(double A,double mass,double mu,double nlq,double nl){
     double lmM=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=A*(lmM)/6.;
     sum[2]=A*A*(7./24. + (19.*lmM)/24. + (lmM*lmM)/36.);
     sum[3]=A*A*A*(58933./124416. + (2.*Zeta2)/3. + (2.*Zeta2*log(2))/9. + 
                    (80507.*Zeta3)/27648. + (8941.*lmM)/1728. +
                    (511.*lmM*lmM)/576. + (lmM*lmM*lmM)/216. +
                    nlq*(-2479./31104. - Zeta2/9. - (409.*lmM)/1728));
     sum[4]=A*A*A*A*
       (39.754401025244204746 + (47039*pow(lmM,2))/3456 + 
	(14149*pow(lmM,3))/10368 + pow(lmM,4)/1296 - 
	(2965705*A4)/54432 + (121  *  A5  )/36 + 
	(697121*Zeta2)/19440 + (49*log(2)*Zeta2)/54 + (243892631*Zeta4)/8709120 - 
	(605*Zeta2*Zeta4)/2688 - (330575*Zeta5)/41472 - (587*Zeta2*log(2))/81 + 
	(2057*Zeta4*log(2))/576 + (2699593*Zeta2*pow(log(2),2))/217728 + 
	(121*Zeta2*pow(log(2),3))/432 + 
	nlq*(-1773073/746496 - (9115*pow(lmM,2))/10368 - 
	     (107*pow(lmM,3))/1728 - (173*A4)/5184 
	     - (557*Zeta2)/162 + (697709*Zeta4)/165888 - 
	     (115*Zeta5)/576 - (22*Zeta2*log(2))/81 + 
	     (1709*Zeta2*pow(log(2),2))/20736 - 
	     (173*pow(log(2),4))/124416 + 
	     lmM*(-1140191/373248 - (47*Zeta2)/54 - 
		  (2*log(2)*Zeta2)/27 - (132283*Zeta3)/82944) 
	     - (4756441*Zeta3)/995328) + 
	nlq*nlq*(140825/1492992 + (493*pow(lmM,2))/20736 + 
		 lmM*(1679/186624 + Zeta2/27) + 
		 (13*Zeta2)/162 + (19*Zeta3)/1728) - (1439*Zeta2*Zeta3)/216 + 
	lmM*(21084715/746496 + (35*Zeta2)/9 + (35*log(2)*Zeta2)/27 + 
	     (2922161*Zeta3)/165888));

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsUpOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsUpOS(double als, double massth, double muth, int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fas6to5os(als/Pi, massth, muth, Nf,nl));
     return als*erg;
}


double CRunDec::fas6to5ms(double A,double mass,double mu,double nlq,double nl){
     double lmm=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[4];
     sum[0]=  1.;
     sum[1]=A*(lmm)/6.;
     sum[2]=A*A*(-11./72. + (11.*lmm)/24. + lmm*lmm/36.);
     sum[3]=A*A*A*(-564731./124416. + (2645.*lmm)/1728. + (167.*lmm*lmm)/576. + lmm*lmm*lmm/216. + 
   					(2633./31104. - (67.*lmm)/576. + lmm*lmm/36.)*nlq + (82043.*Zeta3)/27648.);

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsUpMS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsUpMS(double als, double massth, double muth, int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fas6to5ms(als/Pi, massth, muth, Nf,nl));
     return als*erg;
}


double CRunDec::fas6to5si(double A,double mass,double mu,double nlq, double nl){
     double lmmu=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[4];
     sum[0] = 1.0;
     sum[1] = (A*lmmu)/6.;
     sum[2] = A*A*(-11./72. + (19.*lmmu)/24. + lmmu*lmmu/36.);
     sum[3] = A*A*A*(-564731./124416. + (2191.*lmmu)/576. + (511.*lmmu*lmmu)/576. + lmmu*lmmu*lmmu/216. + 
              (2633./31104. - (281.*lmmu)/1728.)*nlq + (82043.*Zeta3)/27648.);

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsUpSI(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsUpSI(double als, double massth, double muth, int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fas6to5si(als/Pi, massth, muth, Nf,nl));
     return als*erg;
}

double CRunDec::fas5to6si(double A,double mass,double mu,double nlq,double nl){
     double lmmu=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[4];
     sum[0] = 1;
     sum[1] = -(A*lmmu)/6.;
     sum[2] = A*A*(11./72. - (19.*lmmu)/24. + lmmu*lmmu/36.); 
     sum[3] = A*A*A*(564731./124416. - (6793.*lmmu)/1728. - (131.*lmmu*lmmu)/576. - lmmu*lmmu*lmmu/216. + 
              (-2633./31104. + (281.*lmmu)/1728.)*nlq - (82043.*Zeta3)/27648.);

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsDownSI(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsDownSI(double als, double massth, double muth, int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fas5to6si(als/Pi, massth, muth, Nf,nl));
     return als*erg;
}


// Coefficients of eq.(33) of [RunDec]
double CRunDec::fmq6to5os(double A,double mass,double mu,double nlq,double nl){
     double lmM=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=0.;
     sum[2]=-A*A*(89./432. - (5.*lmM)/36. + lmM*lmM/12.);
     sum[3]=-A*A*A*(1871./2916. - (4.*A4)/9. - ln2*ln2*ln2*ln2/54. + (121.*lmM)/2592. + (319.*lmM*lmM)/432. + (29.*lmM*lmM*lmM)/216. + (1327.*nlq)/11664.
			- (53.*lmM*nlq)/432. - (lmM*lmM*lmM*nlq)/108. + (ln2*ln2*Zeta2)/9. - (407.*Zeta3)/864. - (5.*lmM*Zeta3)/6. - (2.*nlq*Zeta3)/27. + (103.*Zeta4)/72.);
     sum[4]=-A*A*A*A*(122466970229./3292047360. - (1924649.*A4)/181440. - (118.*A5)/27. - (1924649.*ln2*ln2*ln2*ln2)/4354560. +
			(59.*ln2*ln2*ln2*ln2*ln2)/1620. - (1126487.*lmM)/373248. - (31.*A4*lmM)/9.-(31.*ln2*ln2*ln2*ln2*lmM)/216. + (104803.*lmM*lmM)/
    10368. + (1403.*lmM*lmM*lmM)/648. + (305.*lmM*lmM*lmM*lmM)/1152. - (2403419.*nlq)/
    746496. + (49.*A4*nlq)/108. + (4.*A5*nlq)/9. + (49.*ln2*ln2*ln2*ln2*nlq)/2592. - (
    ln2*ln2*ln2*ln2*ln2*nlq)/270. + (7045.*lmM*nlq)/31104. + (2.*A4*lmM*nlq)/9.+ 
    1./108.*ln2*ln2*ln2*ln2*lmM*nlq - (9601.*lmM*lmM*nlq)/10368. - (47.*lmM*lmM*lmM*nlq)/
    288. - (5.*lmM*lmM*lmM*lmM*nlq)/144. + (17671.*nlq*nlq)/124416. - (3401.*lmM*nlq*nlq)/
    46656. + (31.*lmM*lmM*nlq*nlq)/1296. + (lmM*lmM*lmM*lmM*nlq*nlq)/864. - (5.*Zeta2)/9. - (
    5.*ln2*Zeta2)/27. + (1924649.*ln2*ln2*Zeta2)/725760. - (59.*ln2*ln2*ln2*Zeta2)/162. + (
    2.*lmM*Zeta2)/3. + (2.*ln2*lmM*Zeta2)/9. + 31./36.*ln2*ln2*lmM*Zeta2 + (5.*nlq*Zeta2)/
    54. - 49./432.*ln2*ln2*nlq*Zeta2 + 1./27.*ln2*ln2*ln2*nlq*Zeta2 - (lmM*nlq*Zeta2)/9. - 
    1./18.*ln2*ln2*lmM*nlq*Zeta2 - (441628331.*Zeta3)/40642560. - (420877.*lmM*Zeta3)/
    27648. - (155.*lmM*lmM*Zeta3)/48. - (1075.*nlq*Zeta3)/1728. + (221.*lmM*nlq*Zeta3)/
    576. - (5.*nlq*nlq*Zeta3)/864. + 7./108.*lmM*nlq*nlq*Zeta3 + (1061.*Zeta3*Zeta3)/576. + (
    16187201.*Zeta4)/580608. - (725.*ln2*Zeta4)/108. + (4123.*lmM*Zeta4)/288. - (
    1225.*nlq*Zeta4)/3456. + (49.*ln2*nlq*Zeta4)/72. - (163.*lmM*nlq*Zeta4)/144. - (
    7.*nlq*nlq*Zeta4)/96. - (59015.*Zeta5)/1728. + (575.*lmM*Zeta5)/72. + (497.*nlq*Zeta5)/
    288. - (3935.*Zeta6)/432.);

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

double CRunDec::fmq6to5ms(double A,double mass,double mu,double nlq,double nl){
     double lmM=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=0.;
     sum[2]=A*A*(-89./432. + (5.*lmM)/36. - (lmM*lmM)/12.);
     sum[3]=A*A*A*(-2951./2916. - Zeta2*log(2)*ln2*ln2/9. + ln2*ln2*ln2*ln2/54. + 407.*Zeta3/864. - 103.*Zeta4/72.
                   + 4.*A4/9. + (311./2592. + 5.*Zeta3/6.)*lmM - 175*lmM*lmM/432.
                   - 29.*lmM*lmM*lmM/216. + nlq*(-1327./11664. + 2.*Zeta3/27.
                   + 53.*lmM/432. + lmM*lmM*lmM/108.));
     sum[4]=A*A*A*A*(-131968227029./3292047360. + 1924649.*ln2*ln2*ln2*ln2/4353450. - 59.*ln2*ln2*ln2*ln2*ln2/1620. - 1924649.*Zeta2*ln2*ln2/725760.
                     + 59.*Zeta2*ln2*ln2*ln2/162. + 353193131.*Zeta3/40642560. - 1061.*Zeta3*Zeta3/576. - 16187201.*Zeta4/580608.
                     + 725.*Zeta4*ln2/108. + 59015.*Zeta5/1728. + 3935.*Zeta6/432. + 1924649.*A4/181440. + 118.*A5/27.
                     + lmM*(2810855./373248. + 31.*ln2*ln2*ln2*ln2/216. - 31.*Zeta2*ln2*ln2/36. + 373261.*Zeta3/27648. - 4123.*Zeta4/288.
                     - 575.*Zeta5/72. + 31.*A4/9.) + lmM*lmM*(-51163./10368. + 155.*Zeta3/48.) - 301.*lmM*lmM*lmM/324. - 305.*lmM*lmM*lmM*lmM/1152.
                     + nlq*(2261435./746496. - 49.*ln2*ln2*ln2*ln2/2592. + ln2*ln2*ln2*ln2*ln2/270. + 49.*Zeta2*ln2*ln2/432. - Zeta2*ln2*ln2*ln2/27.
                     + 1075.*Zeta3/1728. + 1225.*Zeta4/3456. - 49.*Zeta4*ln2/72. - 497.*Zeta5/288. - 49*A4/108. - 4.*A5/9. +lmM*(
                     - 16669./31104. - ln2*ln2*ln2*ln2/108. + Zeta2*ln2*ln2/18. - 221.*Zeta3/576. + 163.*Zeta4/144. - 2.*A4/9.) + 7825.*lmM*lmM/10368.
                     + 23.*lmM*lmM*lmM/288. + 5.*lmM*lmM*lmM*lmM/144.) + nlq*nlq*(-17671./124416 + 5.*Zeta3/864. + 7.*Zeta4/96. 
                     + lmM*(3401./46656. - 7.*Zeta3/108.) - 31.*lmM*lmM/1296. - lmM*lmM*lmM*lmM/864.));
     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsUpOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecMqUpOS(double mq, double als, double massth, double muth, 
                           int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN 
     }
     double erg=(this->fmq6to5os(als/Pi, massth, muth, Nf+1,nl));
     return mq*erg;
}

double CRunDec::DecMqUpMS(double mq, double als, double massth, double muth, 
                           int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN 
     }
     double erg=(this->fmq6to5ms(als/Pi, massth, muth, Nf+1,nl));
     return mq*erg;
}

// Coefficients of eq.(30) of [RunDec]
double CRunDec::fmq5to6os(double A,double mass,double mu,double nlq,double nl){
     double lmM=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=0.;
     sum[2]=A*A*(89./432. - (5.*lmM)/36. + lmM*lmM/12.);
     sum[3]=A*A*A*(1871./2916. - (4.*A4)/9. - ln2*ln2*ln2*ln2/54. + (121.*lmM)/2592. + (319.*lmM*lmM)/432. + (29.*lmM*lmM*lmM)/216. + (1327.*nlq)/11664.
			- (53.*lmM*nlq)/432. - (lmM*lmM*lmM*nlq)/108. + (ln2*ln2*Zeta2)/9. - (407.*Zeta3)/864. - (5.*lmM*Zeta3)/6. - (2.*nlq*Zeta3)/27. + (103.*Zeta4)/72.);
     sum[4]=A*A*A*A*(122466970229./3292047360. - (1924649.*A4)/181440. - (118.*A5)/27. - (1924649.*ln2*ln2*ln2*ln2)/4354560. +
			(59.*ln2*ln2*ln2*ln2*ln2)/1620. - (1126487.*lmM)/373248. - (31.*A4*lmM)/9.-(31.*ln2*ln2*ln2*ln2*lmM)/216. + (104803.*lmM*lmM)/
    10368. + (1403.*lmM*lmM*lmM)/648. + (305.*lmM*lmM*lmM*lmM)/1152. - (2403419.*nlq)/
    746496. + (49.*A4*nlq)/108. + (4.*A5*nlq)/9. + (49.*ln2*ln2*ln2*ln2*nlq)/2592. - (
    ln2*ln2*ln2*ln2*ln2*nlq)/270. + (7045.*lmM*nlq)/31104. + (2.*A4*lmM*nlq)/9.+ 
    1./108.*ln2*ln2*ln2*ln2*lmM*nlq - (9601.*lmM*lmM*nlq)/10368. - (47.*lmM*lmM*lmM*nlq)/
    288. - (5.*lmM*lmM*lmM*lmM*nlq)/144. + (17671.*nlq*nlq)/124416. - (3401.*lmM*nlq*nlq)/
    46656. + (31.*lmM*lmM*nlq*nlq)/1296. + (lmM*lmM*lmM*lmM*nlq*nlq)/864. - (5.*Zeta2)/9. - (
    5.*ln2*Zeta2)/27. + (1924649.*ln2*ln2*Zeta2)/725760. - (59.*ln2*ln2*ln2*Zeta2)/162. + (
    2.*lmM*Zeta2)/3. + (2.*ln2*lmM*Zeta2)/9. + 31./36.*ln2*ln2*lmM*Zeta2 + (5.*nlq*Zeta2)/
    54. - 49./432.*ln2*ln2*nlq*Zeta2 + 1./27.*ln2*ln2*ln2*nlq*Zeta2 - (lmM*nlq*Zeta2)/9. - 
    1./18.*ln2*ln2*lmM*nlq*Zeta2 - (441628331.*Zeta3)/40642560. - (420877.*lmM*Zeta3)/
    27648. - (155.*lmM*lmM*Zeta3)/48. - (1075.*nlq*Zeta3)/1728. + (221.*lmM*nlq*Zeta3)/
    576. - (5.*nlq*nlq*Zeta3)/864. + 7./108.*lmM*nlq*nlq*Zeta3 + (1061.*Zeta3*Zeta3)/576. + (
    16187201.*Zeta4)/580608. - (725.*ln2*Zeta4)/108. + (4123.*lmM*Zeta4)/288. - (
    1225.*nlq*Zeta4)/3456. + (49.*ln2*nlq*Zeta4)/72. - (163.*lmM*nlq*Zeta4)/144. - (
    7.*nlq*nlq*Zeta4)/96. - (59015.*Zeta5)/1728. + (575.*lmM*Zeta5)/72. + (497.*nlq*Zeta5)/
    288. - (3935.*Zeta6)/432.);

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsUpOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecMqDownOS(double mq, double als, double massth, double muth, 
                           int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fmq5to6os(als/Pi, massth, muth, Nf,nl));
     return mq*erg;
}

// Function double CRunDec::AlL2AlH(double asl,double mu1,TriplenfMmu decpar[],
//                          double mu2, int nl)
double CRunDec::AlL2AlH(double asl,double mu1,TriplenfMmu decpar[],double mu2, 
                         int nl){
     int n=0;
     int help;
     double help2;
     double asini=asl;
     double muini=mu1;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf>decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=3-n+2;i<=3;i++){
       if(decpar[i].nf-decpar[i-1].nf!=1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."<<endl;
         RETURN  
       }
     }
     double erg1,erg2;
     int i;
     for(i=3-n+1;i<=3;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf-1,nl);
       erg2= DecAsUpOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       asini=erg2;
       muini=decpar[i].muth; 
     }
     double alpha= AlphasExact(asini,muini,mu2,decpar[i-1].nf,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return alpha;
}

// Function double CRunDec::AlH2AlL(double ash,double mu1,TriplenfMmu decpar[],
//                          double mu2, int nl)
double CRunDec::AlH2AlL(double ash,double mu1,TriplenfMmu decpar[],double mu2, 
                         int nl){
     int n=0;
     int help;
     double help2;
     double asini=ash;
     double muini=mu1;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf<decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=1;i<=n-1;i++){
       if(decpar[i].nf-decpar[i-1].nf!=-1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT.";
         RETURN  
       }
     }
     double erg1,erg2;
     int i;
     for(i=0;i<=n-1;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf,nl);
       erg2= DecAsDownOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       asini=erg2;
       muini=decpar[i].muth; 
     }
     double alpha= AlphasExact(asini,muini,mu2,decpar[i-1].nf-1,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return alpha;
}

// Function double CRunDec::mL2mH(double mql,double asl,double mu1,
//                         TriplenfMmu decpar[], double mu2, int nl)
double CRunDec::mL2mH(double mql,double asl,double mu1,TriplenfMmu decpar[],
                       double mu2, int nl){
     int n=0;
     int help;
     double help2;
     double asini=asl;
     double muini=mu1;
     double mqini=mql;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf>decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=3-n+2;i<=3;i++){
       if(decpar[i].nf-decpar[i-1].nf!=1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT.";
         RETURN  
       }
     }
     
     double erg1,erg2,erg3, erg4, erg5, erg6;
     int i;
     for(i=3-n+1;i<=3;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf-1,nl);
       erg3= mMS2mMS(mqini,asini,erg1,decpar[i].nf-1,nl);
       erg2= DecAsUpOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       erg4=DecMqUpOS(erg3,erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       asini=erg2;
       mqini=erg4;
       muini=decpar[i].muth; 
     }
     erg5= AlphasExact(asini,muini,mu2,decpar[i-1].nf,nl);
     erg6= mMS2mMS(mqini,asini,erg5,decpar[i-1].nf,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return erg6;
}

// Function double CRunDec::mH2mL(double mqh,double ash,double mu1,
//                         TriplenfMmu decpar[], double mu2, int nl)
double CRunDec::mH2mL(double mqh,double ash,double mu1,TriplenfMmu decpar[],
                      double mu2, int nl){
     int n=0;
     int help;
     double help2;
     double asini=ash;
     double muini=mu1;
     double mqini=mqh;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf<decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=1;i<=n-1;i++){
       if(decpar[i].nf-decpar[i-1].nf!=-1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT.";
         RETURN  
       }
     }
     double erg1,erg2,erg3,erg4,erg5,erg6;
     int i;
     for( i=0;i<=n-1;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf,nl);
       erg3= mMS2mMS(mqini,asini,erg1,decpar[i].nf,nl);
       erg2= DecAsDownOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       erg4=DecMqDownOS(erg3,erg1,decpar[i].Mth,decpar[i].muth,
            decpar[i].nf-1,nl);
       asini=erg2;
       mqini=erg4;
       muini=decpar[i].muth; 
     }
     erg5= AlphasExact(asini,muini,mu2,decpar[i-1].nf-1,nl);
     erg6 =mMS2mMS(mqini,asini,erg5,decpar[i-1].nf-1,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return erg6;
}


// In the following the functions above are overloaded w.r.t. to an
// additional argument, n_f (number of active flavours).
// Use SetConstants(nf)
// Functions depending on the 4-loop MS-OS-massrelation are also overloaded w.r.t to the multiplicative error, which is set to 1.0 in this case
double CRunDec::LamExpl(double AlphaS, double Mu, int nf, int nl){
     SetConstants(nf);
     return (this->LamExpl(AlphaS, Mu, nl));      
}
  
double CRunDec::LamImpl(double AlphaS, double Mu,int nf,int nl){
     SetConstants(nf);
     return (this->LamImpl(AlphaS, Mu, nl));
}
  
double CRunDec::AlphasLam(double Lambda, double Mu,int nf, int nl){
     SetConstants(nf);
     return (this->AlphasLam(Lambda, Mu, nl));
}
  
double CRunDec::AlphasExact(double AlphaS0, double Mu0, double MuEnd, 
                                int nf,int nl){
     SetConstants(nf);
     return (this->AlphasExact(AlphaS0,Mu0,MuEnd,nl));
}
  
double CRunDec::mMS2mMS(double Mu0, double AlphaSEnd, double AlphaS0,int nf, 
                           int nl){
     SetConstants(nf);
     return (this->mMS2mMS(Mu0, AlphaSEnd, AlphaS0, nl));
}
  
AsmMS CRunDec::AsmMSrunexact(double mMu, double AlphaS0, double Mu0, 
                              double MuEnd,int nf, int nl){
     SetConstants(nf);
     return (this->AsmMSrunexact(mMu, AlphaS0, Mu0, MuEnd, nl));
}

double CRunDec::mMS2mOS(double MS, double mq[], double asmu,double mu,int nf,
                           int nl, double err){
     SetConstants(nf);
     return (this->mMS2mOS(MS, mq, asmu, mu, nl, err));
} 
  
double CRunDec::mOS2mMS(double mOS, double mq[],double asmu,double Mu,int nf,
                           int nl, double err){
     SetConstants(nf);
     return (this->mOS2mMS(mOS, mq, asmu, Mu, nl, err));
}
    
double CRunDec::mMS2mOS(double MS, double mq[], double asmu,double mu,int nf,
                           int nl){
     SetConstants(nf);
     return (this->mMS2mOS(MS, mq, asmu, mu, nl, 1.0));
} 
  
double CRunDec::mOS2mMS(double mOS, double mq[],double asmu,double Mu,int nf,
                           int nl){
     SetConstants(nf);
     return (this->mOS2mMS(mOS, mq, asmu, Mu, nl, 1.0));
}

double CRunDec::mMS2mOS(double MS, double mq[], double asmu,double mu,
                           int nl){
     return (this->mMS2mOS(MS, mq, asmu, mu, nl, 1.0));
} 
  
double CRunDec::mOS2mMS(double mOS, double mq[],double asmu,double Mu,
                           int nl){
     return (this->mOS2mMS(mOS, mq, asmu, Mu, nl, 1.0));
}
  
double CRunDec::mMS2mSI(double mMS, double asmu, double mu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mSI(mMS, asmu, mu, nl));
}
  
double CRunDec::mRI2mMS(double mRI, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mRI2mMS(mRI, asmu, nl));
}
  
double CRunDec::mMS2mRGI(double mMS, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mRGI(mMS, asmu, nl));
}
  
double CRunDec::mRGI2mMS(double mRGI, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mRGI2mMS(mRGI, asmu, nl));
}
  
double CRunDec::mOS2mSI(double mOS, double mq[], double asM,int nf, int nl, double err){
     SetConstants(nf);
     return (this->mOS2mSI(mOS, mq, asM, nl, err));
}

double CRunDec::mOS2mSI(double mOS, double mq[], double asM,int nf, int nl){
     SetConstants(nf);
     return (this->mOS2mSI(mOS, mq, asM, nl, 1.0));
}
 
double CRunDec::mOS2mSI(double mOS, double mq[], double asM, int nl){
     return (this->mOS2mSI(mOS, mq, asM, nl, 1.0));
}
 
double CRunDec::mOS2mMSrun(double mOS, double mq[], double asmu, double mu,
                              int nf, int nl){
     SetConstants(nf);
     return (this->mOS2mMSrun(mOS, mq, asmu, mu, nl));
}
   
double CRunDec::mMS2mOSrun(double mMS, double mq[], double asmu, double mu,
                              int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mOSrun(mMS, mq, asmu, mu, nl));
}
   
double CRunDec::mMS2mRI(double mMS, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mRI(mMS, asmu, nl));
}
  
double CRunDec::mOS2mMSit(double mOS, double mq[], double asmu, double mu,
                             int nf,int nl){
     SetConstants(nf);
     return (this->mOS2mMSit(mOS, mq, asmu, mu, nl));
}
  
double CRunDec::mMS2mRGImod(double mMS, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mRGImod(mMS, asmu, nl));
}       

double CRunDec::mMS2mPS(double mMS, double mq[], double asmu, double mu, double muf, int nl, int nloops) {
     return (this->mMS2mPS(mMS, mq, asmu, mu, muf, nl, nloops, 1.0));
}

double CRunDec::mPS2mMS(double mPS, double mq[], double asmu, double mu, double muf, int nl, int nloops) {
     return (this->mPS2mMS(mPS, mq, asmu, mu, muf, nl, nloops, 1.0));
}

double CRunDec::mPS2mSI(double mPS, double mq[], double (*as)(double), double muf, int nl, int nloops) {
     return (this->mPS2mSI(mPS, mq, as, muf, nl, nloops, 1.0));
}

double CRunDec::mMS2m1S(double mMS, double mq[], double asmu, double mu, int nl, int nloops) {
     return (this->mMS2m1S(mMS, mq, asmu, mu, nl, nloops, 1.0));
}

double CRunDec::m1S2mMS(double m1S, double mq[], double asmu, double mu, int nl, int nloops) {
     return (this->m1S2mMS(m1S, mq, asmu, mu, nl, nloops, 1.0));
}

double CRunDec::m1S2mSI(double m1S, double mq[], double (*as)(double), int nl, int nloops) {
     return (this->m1S2mSI(m1S, mq, as, nl, nloops, 1.0));
}

double CRunDec::mOS2mRS(double mOS, double mq[], double asmu, double mu, double nuf, int nl, int nloops) {
     return (this->mOS2mRS(mOS, mq, asmu, mu, nuf, nl, nloops, false));
}

double CRunDec::mOS2mRSp(double mOS, double mq[], double asmu, double mu, double nuf, int nl, int nloops) {
     return (this->mOS2mRS(mOS, mq, asmu, mu, nuf, nl, nloops, true));
}

double CRunDec::mMS2mRS(double mMS, double mq[], double asmu, double mu, double nuf, int nl, int nloops, double err) {
     return (this->mMS2mRS(mMS, mq, asmu, mu, nuf, nl, nloops, err, false));
}

double CRunDec::mMS2mRSp(double mMS, double mq[], double asmu, double mu, double nuf, int nl, int nloops, double err) {
     return (this->mMS2mRS(mMS, mq, asmu, mu, nuf, nl, nloops, err, true));
}

double CRunDec::mRS2mMS(double mRS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err) {
     return (this->mRS2mMS(mRS, mq, asmu, mu, muf, nl, nloops, err, false));
}

double CRunDec::mRSp2mMS(double mRS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err) {
     return (this->mRS2mMS(mRS, mq, asmu, mu, muf, nl, nloops, err, true));
}

double CRunDec::mRS2mSI(double mRS, double mq[], double (*as)(double), double muf, int nl, int nloops, double err) {
     return (this->mRS2mSI(mRS, mq, as, muf, nl, nloops, err, false));
}

double CRunDec::mRSp2mSI(double mRS, double mq[], double (*as)(double), double muf, int nl, int nloops, double err) {
     return (this->mRS2mSI(mRS, mq, as, muf, nl, nloops, err, true));
}

double CRunDec::mMS2mRS(double mMS, double mq[], double asmu, double mu, double muf, int nl, int nloops) {
     return (this->mMS2mRS(mMS, mq, asmu, mu, muf, nl, nloops, 1.0));
}

double CRunDec::mRS2mMS(double mRS, double mq[], double asmu, double mu, double muf, int nl, int nloops) {
     return (this->mRS2mMS(mRS, mq, asmu, mu, muf, nl, nloops, 1.0));
}

double CRunDec::mRS2mSI(double mRS, double mq[], double (*as)(double), double muf, int nl, int nloops) {
     return (this->mRS2mSI(mRS, mq, as, muf, nl, nloops, 1.0));
}

double CRunDec::mMS2mRSp(double mMS, double mq[], double asmu, double mu, double muf, int nl, int nloops) {
     return (this->mMS2mRSp(mMS, mq, asmu, mu, muf, nl, nloops, 1.0));
}

double CRunDec::mRSp2mMS(double mRS, double mq[], double asmu, double mu, double muf, int nl, int nloops) {
     return (this->mRSp2mMS(mRS, mq, asmu, mu, muf, nl, nloops, 1.0));
}

double CRunDec::mRSp2mSI(double mRS, double mq[], double (*as)(double), double muf, int nl, int nloops) {
     return (this->mRSp2mSI(mRS, mq, as, muf, nl, nloops, 1.0));
}

double CRunDec::DecAsDownOS(double als, double massth, double muth,int nf, 
                               int nl){
     SetConstants(nf);
     return (this->DecAsDownOS(als, massth, muth, nl));
}
  
double CRunDec::DecAsUpOS(double als, double massth, double muth, int nf,
                             int nl){
     SetConstants(nf);
     return (this->DecAsUpOS(als, massth, muth, nl)); 
}

double CRunDec::DecAsUpMS(double als, double massth, double muth, int nf,
                             int nl){
     SetConstants(nf);
     return (this->DecAsUpMS(als, massth, muth, nl)); 
}

double CRunDec::DecAsUpSI(double als, double massth, double muth, int nf,
                             int nl){
     SetConstants(nf);
     return (this->DecAsUpSI(als, massth, muth, nl)); 
}

double CRunDec::DecAsDownSI(double als, double massth, double muth, int nf,
                             int nl){
     SetConstants(nf);
     return (this->DecAsDownSI(als, massth, muth, nl)); 
}
  
double CRunDec::DecMqUpOS(double mq, double als, double massth, double muth,
                             int nf, int nl){
     SetConstants(nf);
     return (this->DecMqUpOS(mq, als, massth, muth, nl));
}
  
double CRunDec::DecMqDownOS(double mq, double als, double massth, 
                               double muth,int nf, int nl){
     SetConstants(nf);
     return (this->DecMqDownOS(mq, als, massth, muth, nl));
}

