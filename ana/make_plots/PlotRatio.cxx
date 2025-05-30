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
void drawplots(string file1, string file2, string hist11name, string hist12name,string hist21name, string hist22name, bool isNSF)
{

  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  //gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);

  TFile f1(Form("%s",file1.c_str()));
  TFile f2(Form("%s",file2.c_str()));
/*  MnvH1D* num1=(MnvH1D*)f1.Get(hist11name.c_str());
  MnvH1D* den1=(MnvH1D*)f1.Get(hist12name.c_str());
  MnvH1D* num2=(MnvH1D*)f2.Get(hist21name.c_str());
  MnvH1D* den2=(MnvH1D*)f2.Get(hist22name.c_str());
*/
  MnvH2D* num1=(MnvH2D*)f1.Get(hist11name.c_str());
  MnvH2D* den1=(MnvH2D*)f1.Get(hist12name.c_str());
  MnvH2D* num2=(MnvH2D*)f2.Get(hist21name.c_str());
  MnvH2D* den2=(MnvH2D*)f2.Get(hist22name.c_str());


/*
  TH1* num1=h1->ProjectionX(); 
  TH1* den1=h2->ProjectionX(); 
  TH1* num2=h3->ProjectionX(); 
  TH1* den2=h4->ProjectionX(); 
  
  TH1* num1=h1->ProjectionY(); 
  TH1* den1=h2->ProjectionY(); 
  TH1* num2=h3->ProjectionY(); 
  TH1* den2=h4->ProjectionY();*/ 
  // qemc->Scale(madpot/qemcpot ,"width");
 
  num1->Scale(1,"width");
  den1->Scale(1,"width");
  num2->Scale(1,"width");
  num2->Scale(0.421);
  den2->Scale(1,"width");

 

  num1->SetMarkerStyle(kFullCircle);
  num1->SetMarkerSize(0.9);
  num1->SetLineColor(kBlue);
  num1->SetMarkerColor(kBlack);
  num2->SetMarkerStyle(kFullCircle);
  num2->SetMarkerSize(0.9);
  num2->SetLineColor(kBlue);
  num2->SetMarkerColor(kBlack);
  applyStyle((TH1D*)num1);
  applyStyle((TH1D*)den1);
  applyStyle((TH1D*)num2);
  applyStyle((TH1D*)den2);

  den1->SetLineColor(kBlack);
  den1->SetMarkerColor(kBlack);
  den2->SetLineColor(kBlack);
  den2->SetMarkerColor(kBlack);


  TCanvas *c= new TCanvas;
  c->cd();
  
  TLegend* leg=new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(num1, "Numerator", "l");
  leg->AddEntry(den1, "Denominator", "l");

  //num1->SetAxisRange(0,4000,"y");
  //num1->SetAxisRange(0,50,"x");
  num1->SetTitle(Form("%s  -  Iron of all Target ",isNSF? "":""));
  //num->SetTitle(Form("%s  -  Lead of Target 1  -  minervame1L",isNSF? "NSFNuke":"NukeCC"));
  //num->GetXaxis()->SetTitle("True Hadron Energy (GeV)");
  //num1->GetXaxis()->SetTitle("True y");
  //num1->GetXaxis()->SetTitle("Muon Energy (GeV)");
  //num1->GetYaxis()->SetTitle("Events/GeV");

  //num1->Draw("hist");
  //den1->Draw("histsame");
  num2->Draw("histsame");
  leg->Draw("same");

  c->Print(Form("EfficiencyNumDen_%s.png",isNSF? "NSFNuke":"NukeCC"));
  //c->Print(Form("EfficiencyNumDen_%s.C",isNSF? "NSFNuke":"NukeCC"));
 //TFile *f = new TFile("test.root","RECREATE");


  
  //num1->Divide(num,den);
  num2->AddMissingErrorBandsAndFillWithCV(*num1);
  num1->Divide(num1, num2);
  //num1->GetYaxis()->SetTitle("Ratio MAD/NukeCC efficiency script");
  //num1->GetYaxis()->SetTitle("Ratio MasterAnaDev/NukeCC Migration script");
  //num1->GetYaxis()->SetTitle("Data / Mnv_GENIE");
  num1->GetYaxis()->SetTitle("Hadronic Energy (GeV)");
  num1->GetXaxis()->SetTitle("Muon Energy (GeV)");
  num1->SetAxisRange(0.9,1.1,"y");
  num1->SetLineColor(kBlack);
  num1->Draw("colz");
  TFile *f = new TFile("test.root","RECREATE");
  //num1->Draw(" ");
  num1->Write();
  f->Write();
  f->Close();
 /* 
  den1->Divide(den1,den2);
  den1->GetYaxis()->SetTitle("Ratio MAD/NukeCC Efficiency script");
  den1->SetAxisRange(0.5,1.5,"y");
  den1->SetLineColor(kBlack);
  den1->GetXaxis()->SetTitle("Muon Energy (GeV)");
  den1->Draw();*/
  c->Print(Form("Efficiency_%s.png",isNSF? "NSFNuke":"NukeCC"));
  c->Print(Form("Efficiency_%s.C",isNSF? "NSFNuke":"NukeCC"));
  c->Print(Form("Efficiency_%s.root",isNSF? "NSFNuke":"NukeCC"));
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
  bool isNSF=atoi(argv[7]);
  drawplots(argv[1],argv[2],argv[3],argv[4], argv[5], argv[6],isNSF);
  return 0;
}
