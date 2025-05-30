#include "myPlotStyle.h"
//#include "GridCanvas.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
// #include "../util/plot/plot.h"                                                                        
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TEnv.h"
#include "TChain.h"
//#include "TF2.h"
#include "Math/DistFunc.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TFile.h"
#include "localColor.h"
#include "Cintex/Cintex.h"
#include "TVector2.h"

#include "myPlotStyle.h"

//#include "plot.h"
using namespace std;
using namespace PlotUtils;
//void drawplots(string location, string playlist)
void drawplots(string file,string hist1name, string hist2name, bool isNSF)
{

  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  //gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);

  TFile f1(Form("%s",file.c_str()));

  MnvH1D* num=(MnvH1D*)f1.Get(hist1name.c_str());
  MnvH1D* den=(MnvH1D*)f1.Get(hist2name.c_str());

  // qemc->Scale(madpot/qemcpot ,"width");
 
  num->Scale(1,"width");
  den->Scale(1,"width");

 

  num->SetMarkerStyle(kFullCircle);
  num->SetMarkerSize(0.9);
  num->SetLineColor(kBlue);
  num->SetMarkerColor(kBlack);
  applyStyle((TH1D*)num);
  applyStyle((TH1D*)den);

  den->SetLineColor(kBlack);
  den->SetMarkerColor(kBlack);

  TCanvas *c= new TCanvas;
  c->cd();
  
  TLegend* leg=new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(num, "Numerator", "l");
  leg->AddEntry(den, "Denominator", "l");

  num->SetAxisRange(0,4000,"y");
  num->SetAxisRange(0,50,"x");
  num->SetTitle(Form("%s  -  Lead of Target 1  -  minervame1A",isNSF? "NSFNuke":"NukeCC"));
  //num->SetTitle(Form("%s  -  Lead of Target 1  -  minervame1L",isNSF? "NSFNuke":"NukeCC"));
  //num->GetXaxis()->SetTitle("True Hadron Energy (GeV)");
  num->GetXaxis()->SetTitle("True y");
  num->GetYaxis()->SetTitle("Events/GeV");

  num->Draw("hist");
  den->Draw("histsame");
  leg->Draw("same");

  c->Print(Form("EfficiencyNumDen_%s.png",isNSF? "NSFNuke":"NukeCC"));
  //c->Print(Form("EfficiencyNumDen_%s.C",isNSF? "NSFNuke":"NukeCC"));
  
  num->Divide(num,den);
  num->GetYaxis()->SetTitle("Efficiency");
  num->SetAxisRange(0,1.0,"y");

  num->SetLineColor(kBlack);
  num->Draw();
  c->Print(Form("Efficiency_%s.png",isNSF? "NSFNuke":"NukeCC"));
  //c->Print(Form("Efficiency_%s.C",isNSF? "NSFNuke":"NukeCC"));
}


int main(int argc, char* argv[]){

  if(argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./PlotEffs filepath/name numeratorhist denominatorhist isNSF"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    return 0;
  }
  bool isNSF=atoi(argv[4]);
  drawplots(argv[1],argv[2],argv[3],isNSF);
  return 0;
}
