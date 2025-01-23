/*
  CRunDec.h

  Header file for CRunDec.cpp

  Author: Barbara Schmidt (Jan 2012)
          Florian Herren and Matthias Steinhauser (Jan 2016)

*/

/*
  Minor update: (Sep 2016)
  Remove "static" from definition and initialization of constants for Runge-Kutta procedure.
  The C++11 standard supports the new version, so a C++11 compliant compiler and maybe
  the switch "-std=c++11" should be used.
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


#ifndef CRUNDEC_H_
#define CRUNDEC_H_

// Default return statement:
#define RETURN return(0);
// The following might be useful for Windows:
//#define RETURN   system("PAUSE"); exit(1);


// Numerical values for input parameters:
#define asMz 0.1185
#define Mz   91.1876
#define Mt   173.2
#define Mb   4.8
#define Mc   1.5
#define muc  1.279
#define mub  4.163
#define Mtau 1.777

// Some constants:
#define cf 4./3.
#define ca 3.
#define tr 1./2.
#define B4 -1.762800087073770864061897634679818807215137274389016762629478603776
#define A4 0.5174790616738993863307581618988629456223774751413792582443193479770
#define A5 0.5084005792422687074591088492585899413195411256648216487244977963526
#define Pi M_PI
#define Zeta2 (Pi*Pi)/6.
#define Zeta3 1.20205690315959428539973816151144999076498629234049888179227155534
#define Zeta4 (Pi*Pi*Pi*Pi)/90.
#define Zeta5 1.03692775514336992633136548645703416805708091950191281197419267790
#define Zeta6 (Pi*Pi*Pi*Pi*Pi*Pi)/945.
#define Zeta7 1.00834927738192282683979754984979675959986356056523870641728313657160147831735573534609696891385132
#define ln2 log(2)


// Struct for triple {nf, Mth, Muth}:
struct TriplenfMmu{
      int nf;
      double Mth;
      double muth;
};

struct AsmMS{
      double Asexact;
      double mMSexact;
};

// class declaration of CRunDec:
class CRunDec
{
  private:
  // Aux. constants for implicit Runge-Kutta-Procedure:
  const double a2=0.2, a3=0.3, a4=0.6, a5=1., a6=0.875;

  const double b21=0.2, b31=3./40., b32=9./40., b41=0.3, b42=-0.9,
                      b43=6./5.;
  const double b51=-11./54., b52=2.5, b53=-70./27., b54=35./27.;
  const double b61=1631./55296., b62=175./512., b63=575./13824.;
  const double b64=44275./110592., b65=253./4096.;

  const double c1=37./378., c2=0., c3=250./621., c4=125./594., c5=0.;
  const double c6= 512./1771.;

  const double dc1=37./378.-2825./27648., dc2=0.-0.,
                      dc3=250./621.-18575./48384.;
  const double dc4=125./594.-13525./55296., dc5=0.-277./14336.,
                      dc6=512./1771.-0.25;

  // Coefficients for diff. equations:
  double Beta[5], B[5], Gamma[5], C[5], Nf;

  // Define constants (if not already done with constructor):
  void SetConstants(int n);

  // R.h.s. of diff. equations:
  friend double fSetdydx(CRunDec S, double A,int nl);

  friend double fSetdydxa1(CRunDec S,double x, double A);
  friend double fSetdydxM1(CRunDec S,double A, double M);

  friend double fSetdydxa2(CRunDec S,double x, double A);
  friend double fSetdydxM2(CRunDec S,double A, double M);

  friend double fSetdydxa3(CRunDec S,double x, double A);
  friend double fSetdydxM3(CRunDec S,double A, double M);

  friend double fSetdydxa4(CRunDec S,double x, double A);
  friend double fSetdydxM4(CRunDec S,double A, double M);

  friend double fSetdydxa5(CRunDec S,double x, double A);
  friend double fSetdydxM5(CRunDec S,double A, double M);

  // Additional aux. functions:
  int Abbruch(void);
  double fSetAsL(double Lambda, double Mu, int nl, double AlphaS);
  double fSetcx(double x, int nl);
  double fOsFromMs1(double mu, double M);
  double fOsFromMs2(double mu, double M, double nl);
  double fOsFromMs3(double mu, double M, double nl);
  double fOsFromMs4(double mu, double M, double nl, double err);
  double fMsFromOs1(double mu, double M);
  double fMsFromOs2(double mu, double M, double nl);
  double fMsFromOs3(double mu, double M, double nl);
  double fMsFromOs4(double mu, double M, double nl, double err);
  double fZmM(double n);
  double fZmInvM(double n);
  double fDelta(double mOS,double mq[]);
  double fMsFromRi1(void);
  double fMsFromRi2(void);
  double fMsFromRi3(void);
  double fMumFromOs1(void);
  double fMumFromOs2(void);
  double fMumFromOs3(void);
  double fMumFromOs4(double err);
  double fRiFromMs(double alpha, double nl);
  double fMsFromRi(double alpha, double nl);
  double fHelpmOS2mMSit(double MS,double mOS, double mq[], double asmu,
                        double mu, int nl);
  double fas5to6os(double alpha,double mass, double mu,double nlq, double nl);
  double fas6to5os(double alpha,double mass, double mu,double nlq, double nl);
  double fas6to5ms(double alpha,double mass, double mu,double nlq, double nl);
  double fas6to5si(double alpha,double mass, double mu,double nlq, double nl);
  double fas5to6si(double alpha,double mass, double mu,double nlq, double nl);
  double fmq5to6os(double A,double mass,double mu,double nlq,double nl);
  double fmq6to5os(double A,double mass,double mu,double nlq,double nl);
  double fmq6to5ms(double A,double mass,double mu,double nlq,double nl);

  double PSdelta(double asmu, double muf, double mu, int nl, int nloops);
  double E1p(double mOS, double asmu, double mu, int nl, int nloops);
  double exOS2RS(double api, double mmu, double nnuf, int nnl, int nloops);
  double exOS2RSp(double api, double mmu, double nnuf, int nnl, int nloops);
  double mMS2mOSmod(double MS, double mq[4], double asmu, double mu, int nf, int nloops, double err);

  double fRungeKuttaImpl(double &x, double y,double &htry, int nl,
                      double (*f)(CRunDec,double, int));
  double fRKSchritt(double x,double y,double h,double &yerr,
                    double (*f)(CRunDec, double ,double));

  public:
  // constructor:
  CRunDec();
  CRunDec(int);

  // Arrays and structs to store data:
  double mq[4];
  TriplenfMmu nfMmu[4];
  AsmMS AM;

  // Function to obtain current number of active flavours:
  int GetNf();
  // Function to set number of active flavours:
  void SetNf(int nf);

  // Functions for the running of alpha_s:
  double LamExpl(double asmu, double mu, int nloops);
  double LamImpl(double asmu, double mu,int nloops);
  double AlphasLam(double Lambda, double mu, int nloops);
  double AlphasExact(double asmu0, double mu0, double mu1, int nloops);

  // Funktions for the runnung of mq and the
  // various mass definitions:
  double mMS2mMS(double mu0, double asmu0, double asmu1, int nloops);
  double mMS2mOS(double MS, double mq[4], double asmu, double mu, int nloops, double err);
  double mOS2mMS(double mOS, double mq[], double asmu, double mu, int nloops, double err);
  double mMS2mSI(double mMS, double asmu, double mu, int nloops);
  double mRI2mMS(double mRI, double asmu, int nloops);
  double mMS2mRGI(double mMS, double asmu, int nloops);
  double mRGI2mMS(double mRGI, double asmu, int nloops);
  double mOS2mSI(double mOS, double mq[], double asM, int nloops, double err);
  double mOS2mMSrun(double mOS, double mq[], double asmu, double mu, int nloops);
  double mMS2mOSrun(double mMS, double mq[], double asmu, double mu, int nloops);
  double mMS2mRI(double mMS, double asmu, int nloops);
  double mOS2mMSit(double mOS, double mq[], double asmu, double mu, int nloops);
  double mMS2mRGImod(double mMS, double asmu, int nloops);


  double mOS2mPS(double mOS, double mq[], double asmu, double mu, double muf, int nl, int nloops);
  double mMS2mPS(double mMS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err);
  double mPS2mMS(double mPS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err);
  double mPS2mSI(double mPS, double mq[], double (*as)(double), double muf, int nl, int nloops, double err);

  double mOS2m1S(double mOS, double mq[], double asmu, double mu, int nl, int nloops);
  double mMS2m1S(double mMS, double mq[], double asmu, double mu, int nl, int nloops, double err);
  double m1S2mMS(double m1S, double mq[], double asmu, double mu, int nl, int nloops, double err);
  double m1S2mSI(double m1S, double mq[], double (*as)(double), int nl, int nloops, double err);

  double mOS2mRS(double mOS, double mq[], double asmu, double mu, double nuf, int nl, int nloops, bool prime);
  double mMS2mRS(double mMS, double mq[], double asmu, double mu, double nuf, int nl, int nloops, double err, bool prime);
  double mRS2mMS(double mRS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err, bool prime);
  double mRS2mSI(double mRS, double mq[], double (*as)(double), double muf, int nl, int nloops, double err, bool prime);

  // Solve coupled differential equations for alpha_s and mq:
  AsmMS AsmMSrunexact(double mmu, double asmu0, double mu0, double mu1,
                      int nloops);

  // Decoupling relations:
  double DecAsDownOS(double asmu, double massth, double muth, int nloops);
  double DecAsUpOS(double asmu, double massth, double muth, int nloops);
  double DecAsUpMS(double asmu, double massth, double muth, int nloops);
  double DecAsUpSI(double asmu, double massth, double muth, int nloops);
  double DecAsDownSI(double asmu, double massth, double muth, int nloops);
  double DecMqUpOS(double mq, double asmu, double massth, double muth, int nloops);
  double DecMqDownOS(double mq, double asmu, double massth, double muth, int nloops);
  double DecMqUpMS(double mq, double asmu, double massth, double muth, int nloops);
  // Running and decoupling:
  double AlL2AlH(double asl,double mu1,TriplenfMmu decpar[],double mu2, int nloops);
  double AlH2AlL(double ash,double mu1,TriplenfMmu decpar[],double mu2, int nloops);
  double mL2mH(double mql,double asl,double mu1,TriplenfMmu decpar[],double mu2,
               int nloops);
  double mH2mL(double mqh,double ash,double mu1,TriplenfMmu decpar[],double mu2,
               int nloops);

  // Overload functions:
  double LamExpl(double asmu, double mu, int nf, int nloops);
  double LamImpl(double asmu, double mu,int nf,int nloops);
  double AlphasLam(double Lambda, double mu,int nf, int nloops);
  double AlphasExact(double asmu0, double mu0, double mu1, int nf,int nloops);
  double mMS2mMS(double mu0, double asmu1, double asmu0,int nf, int nloops);
  AsmMS AsmMSrunexact(double mmu, double asmu0, double mu0, double mu1,
                       int nf, int nloops);
  double mMS2mOS(double MS, double mq[], double asmu, double mu,int nf, int nloops, double err);
  double mOS2mMS(double mOS, double mq[], double asmu, double mu,int nf,int nloops, double err);
  double mMS2mOS(double MS, double mq[], double asmu, double mu,int nf, int nloops);
  double mOS2mMS(double mOS, double mq[], double asmu, double mu,int nf,int nloops);
  double mMS2mOS(double MS, double mq[4], double asmu, double mu, int nloops);
  double mOS2mMS(double mOS, double mq[], double asmu, double mu, int nloops);
  double mMS2mSI(double mMS, double asmu, double mu,int nf, int nloops);
  double mRI2mMS(double mRI, double asmu,int nf, int nloops);
  double mMS2mRGI(double mMS, double asmu,int nf, int nloops);
  double mRGI2mMS(double mRGI, double asmu,int nf, int nloops);
  double mOS2mSI(double mOS, double mq[], double asM,int nf, int nloops, double err);
  double mOS2mSI(double mOS, double mq[], double asM,int nf, int nloops);
  double mOS2mSI(double mOS, double mq[], double asM, int nloops);
  double mOS2mMSrun(double mOS, double mq[], double asmu, double mu,int nf,
                    int nloops);
  double mMS2mOSrun(double mMS, double mq[], double asmu, double mu,int nf,
                    int nloops);
  double mMS2mRI(double mMS, double asmu,int nf, int nloops);
  double mOS2mMSit(double mOS, double mq[],double asmu,double mu,int nf,int nloops);
  double mMS2mRGImod(double mMS, double asmu,int nf, int nloops);
  double mMS2mPS(double mMS, double mq[], double asmu, double mu, double muf, int nl, int nloops);
  double mPS2mMS(double mPS, double mq[], double asmu, double mu, double muf, int nl, int nloops);
  double mPS2mSI(double mPS, double mq[], double (*as)(double), double muf, int nl, int nloops);
  double mMS2m1S(double mMS, double mq[], double asmu, double mu, int nl, int nloops);
  double m1S2mMS(double m1S, double mq[], double asmu, double mu, int nl, int nloops);
  double m1S2mSI(double m1S, double mq[], double (*as)(double), int nl, int nloops);
  double mOS2mRS(double mOS, double mq[], double asmu, double mu, double nuf, int nl, int nloops);
  double mMS2mRS(double mMS, double mq[], double asmu, double mu, double nuf, int nl, int nloops, double err);
  double mRS2mMS(double mRS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err);
  double mRS2mSI(double mRS, double mq[], double (*as)(double), double muf, int nl, int nloops, double err);
  double mMS2mRS(double mMS, double mq[], double asmu, double mu, double muf, int nl, int nloops);
  double mRS2mMS(double mRS, double mq[], double asmu, double mu, double muf, int nl, int nloops);
  double mRS2mSI(double mRS, double mq[], double (*as)(double), double muf, int nl, int nloops);
  double mOS2mRSp(double mOS, double mq[], double asmu, double mu, double nuf, int nl, int nloops);
  double mMS2mRSp(double mMS, double mq[], double asmu, double mu, double nuf, int nl, int nloops, double err);
  double mRSp2mMS(double mRS, double mq[], double asmu, double mu, double muf, int nl, int nloops, double err);
  double mRSp2mSI(double mRS, double mq[], double (*as)(double), double muf, int nl, int nloops, double err);
  double mMS2mRSp(double mMS, double mq[], double asmu, double mu, double muf, int nl, int nloops);
  double mRSp2mMS(double mRS, double mq[], double asmu, double mu, double muf, int nl, int nloops);
  double mRSp2mSI(double mRS, double mq[], double (*as)(double), double muf, int nl, int nloops);
  double DecAsDownOS(double asmu, double massth, double muth,int nf, int nloops);
  double DecAsUpOS(double asmu, double massth, double muth, int nf, int nloops);
  double DecAsUpMS(double asmu, double massth, double muth, int nf, int nloops);
  double DecAsUpSI(double asmu, double massth, double muth, int nf, int nloops);
  double DecAsDownSI(double asmu, double massth, double muth, int nf, int nloops);
  double DecMqUpOS(double mq, double asmu, double massth, double muth, int nf,
                   int nloops);
  double DecMqDownOS(double mq, double asmu, double massth, double muth,int nf,
                     int nloops);

};

#endif
