#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <functional>

#include "CRunDec.h"

#include "LHAPDF/LHAPDF.h"
#include "/home/tom/Documents/software/software/REvolver-master/code/include/REvolver.h"

int main() {
  // Running of the mass
  // extract alphas
  LHAPDF::PDF* pdf = LHAPDF::mkPDF("NNPDF31_nnlo_as_0118", 0);
  const double MZ = 91.1876;
  const double mb = 4.18;
  const double mt = 162.7;
  const double MH = 125.;
  const double alphas = pdf->alphasQ(MZ);
  std::cout << "alphas(MZ) = " << alphas << std::endl;

  revo::RunPar alphasPar = {5, alphas, MZ};
  revo::RunParV mPar;
  mPar.push_back({5, mb, mb});
  mPar.push_back({6, mt, mt});
  revo::RunParV mPar_empty;
  mPar_empty.push_back({6, mt, mt});

  revo::Core core_test(6, alphasPar, mPar_empty, revo::doubleV({1.}), 4);
  revo::Core core1(6, alphasPar, mPar, revo::doubleV({0.2, 0.1}), 1);
  revo::Core core2(6, alphasPar, mPar, revo::doubleV({0.2, 0.1}), 2);
  revo::Core core3(6, alphasPar, mPar, revo::doubleV({0.2, 0.1}), 3);
  revo::Core core4(6, alphasPar, mPar, revo::doubleV({0.2, 0.1}), 4);
  revo::Core core[4] = {core1, core2, core3, core4};

  std::cout << "alphas(MZ) Revolver = " << core4.alpha()(MZ) << "\n" << std::endl;
  std::cout << "mt(mt) = " << mt << std::endl;
  std::cout << "mt(mt) Revolver = " << core4.masses().mMS(6, mt) << "\n" << std::endl;
  std::cout << "mb(mb) = " << mb << std::endl;
  std::cout << "mb(mb) Revolver = " << core4.masses().mMS(5, mb) << "\n" << std::endl;

  // Pre-defined nf=5:
  CRunDec* pObjnf5 = new CRunDec(5);
  double alphas_secdec = pObjnf5 -> AlphasExact(alphas,MZ,MZ,5,4);
  int nf1 = pObjnf5 -> GetNf();
  std::cout << "alpha_s(MZ) SecDec = " << alphas_secdec << std::endl;

  std::cout << "------------------------------------------------------------" << std::endl;

  // ------------------------------------------------------------

  // Compute mb(mb) from mb(10 GeV)=3.610 GeV and alpha_s(Mz)=0.1189

  CRunDec cmpasmb(6);
  double mu,asmb, asmt;

  asmb = cmpasmb.AlphasExact(alphas,MZ,mb,5,4);
  asmt = cmpasmb.AlphasExact(alphas,MZ,mt,5,4);
  AsmMS asmMS = cmpasmb.AsmMSrunexact(mb, asmb, mb, MH/2., 6, 4);

  std::cout << "mb(MH/2) SecDec = " << asmMS.mMSexact << std::endl;
  std::cout << "mb(MH/2) REvolver = " << core4.masses().mMS(5, MH/2.) << std::endl;

  std::cout << "------------------------------------------------------------" << std::endl;

  std::function<double(double)> running_mass_RunDec = [ &cmpasmb, &asmt, &mt](double muR) {
    if(abs(muR - mt) > 15) {
      double mt_new = cmpasmb.AsmMSrunexact(mt, asmt, mt, muR, 5, 4).mMSexact;
      return mt_new;
    }
    else {
      double mt1 = cmpasmb.AsmMSrunexact(mt, asmt, mt, mt - 15., 5, 4).mMSexact;
      double mt2 = cmpasmb.AsmMSrunexact(mt, asmt, mt, mt - 14., 5, 4).mMSexact;
      double mt3 = cmpasmb.AsmMSrunexact(mt, asmt, mt, mt + 14., 5, 4).mMSexact;
      double mt4 = cmpasmb.AsmMSrunexact(mt, asmt, mt, mt + 15., 5, 4).mMSexact;

      double y1 = mt1;
      double y2 = mt4;
      double x1 = mt - 15.;
      double x2 = mt + 15.;
      double dy1 = (mt2 - mt1)/1.;
      double dy2 = (mt4 - mt3)/1.;
      double a0 = -((dy2*std::pow(x1,3)*x2 + dy1*std::pow(x1,2)*std::pow(x2,2) -
                     dy2*std::pow(x1,2)*std::pow(x2,2) - dy1*x1*std::pow(x2,3) -
                     3*x1*std::pow(x2,2)*y1 + std::pow(x2,3)*y1 - std::pow(x1,3)*y2 +
                     3*std::pow(x1,2)*x2*y2)/std::pow(x1 - x2,3));
      double a1 = -((-(dy2*std::pow(x1,3)) - 2*dy1*std::pow(x1,2)*x2 -
                       dy2*std::pow(x1,2)*x2 + dy1*x1*std::pow(x2,2) +
                       2*dy2*x1*std::pow(x2,2) + dy1*std::pow(x2,3) + 6*x1*x2*y1 -
                       6*x1*x2*y2)/std::pow(x1 - x2,3));
      double a2 = -((dy1*std::pow(x1,2) + 2*dy2*std::pow(x1,2) + dy1*x1*x2 - dy2*x1*x2 -
                     2*dy1*std::pow(x2,2) - dy2*std::pow(x2,2) - 3*x1*y1 - 3*x2*y1 +
                     3*x1*y2 + 3*x2*y2)/std::pow(x1 - x2,3));
      double a3 = -((-(dy1*x1) - dy2*x1 + dy1*x2 + dy2*x2 + 2*y1 - 2*y2)/
                     std::pow(x1 - x2,3));
      return a0 + a1*muR + a2*std::pow(muR, 2) + a3*std::pow(muR, 3);
    }

  };

  std::vector<double> muRs[4], alphas_running[4], mt_running[4], mb_running[4], alphas_running_LHAPDF[4], mb_running_RunDec, mt_running_RunDec, alphas_running_RunDec;

  for(int i = 0; i < 4; i++) {
    double muR = 5.;
    while(muR < 1000) {
      muRs[i].push_back(muR);
      alphas_running[i].push_back(core[i].alpha()(muR));
      alphas_running_LHAPDF[i].push_back(pdf->alphasQ(muR));
      //alphas_running[i].push_back(core_test.alpha()(muR) - (pdf->alphasQ(muR)));
      mt_running[i].push_back(core[i].masses().mMS(6, muR));
      mb_running[i].push_back(core[i].masses().mMS(5, muR));
      if(i == 0){
        mb_running_RunDec.push_back(cmpasmb.AsmMSrunexact(mb, asmb, mb, muR, 5, 4).mMSexact);
        mt_running_RunDec.push_back(running_mass_RunDec(muR));
        //mt_running_RunDec.push_back(cmpasmb.AsmMSrunexact(mt, asmt, mt, muR, 5, 4).mMSexact);
        alphas_running_RunDec.push_back(cmpasmb.AlphasExact(alphas,MZ,muR,5,4));
      }
      //mb_running[i].push_back(0);

      muR += .1;
    }
  }


  for(int l = 0; l < 4; l++) {
    std::ofstream outfile;
    outfile.open ("running" + std::to_string(l + 1) + "l.txt");
    for(int i = 0; i < alphas_running[l].size(); i++) {
      //std::cout << muRs[i] << "\t" << alphas_running[i] << "\t" << mt_running[i] << "\t" << mb_running[i] << std::endl;
      outfile << muRs[l][i] << "," << alphas_running[l][i] << "," << mt_running[l][i] << "," << mb_running[l][i] << ","
              << alphas_running_LHAPDF[l][i] << "," << mb_running_RunDec[i] << "," << mt_running_RunDec[i] << "," << alphas_running_RunDec[i] << std::endl;
    }
    outfile.close();
  }

}