#include <bits/stdc++.h>
#include "TClass.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TKey.h"

#include <iostream>
#include <vector>
#include <math.h>
#include <iostream>

#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TEnv.h"
#include "TChain.h"
#include "TF2.h"
#include "Math/DistFunc.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TFile.h"
//#include "localColor.h"
#include "PlotUtils/MnvColors.h"
#include "Cintex/Cintex.h"
#include "plot.h"

#include "myPlotStyle.h"


//////////
// Adapted from code originally by Deepika Jena
// Reintroduced dependency of plot.h as these methods are not specific to one macro and will be used by at minimum 3 for the tki-ptmu analysis
//////////

using namespace std;
#include "vartex.h"
std::map<string, string> varname;

using namespace PlotUtils;


//_____________________________________________________________
void makePlots(bool doMultipliers,bool doGenies,bool doRatio,string location,string varx, string vary, string iiter)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();

  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  //  gStyle->SetErrorX(0);Disabled because we need the widths to show up because of the non-grid nature
  gStyle->SetEndErrorSize(2);

  string var = varx+"_"+vary;
  string varspecfilename = var;
  string varspec = var;
  if(vary=="ptmu") varspecfilename=varx;
  cout << varspecfilename << endl;
  string filename = Form("%s_CV/CrossSection_iter_%s_CombinedPlaylists.root_%s.root",location.c_str(),iiter.c_str(), varspecfilename.c_str());//Final result
  string filenamenolowrec = Form("%s_pion_rpa/CrossSection_iter_%s_CombinedPlaylists.root_%s.root",location.c_str(),iiter.c_str(), varspecfilename.c_str());//Final result
  cout << filename << endl;


  
  TFile f1(filename.c_str());
  TFile f2(filenamenolowrec.c_str());
  //f1.ls();

  MnvH1D* dataMnv           = (MnvH1D*)f1.Get(Form("h_%s_data_nobck_unfold_effcor_cross_section",varspec.c_str()));
  MnvH1D* mcMnv             = (MnvH1D*)f1.Get(Form("h_%s_mc_nobck_unfold_effcor_cross_section",varspec.c_str()));
  MnvH1D* nomGenieMnv       = (MnvH1D*)f1.Get(Form("h_%s_mc_nobck_unfold_effcor_cross_section",varspec.c_str()));//2
  MnvH1D* bestGenieMnv      = (MnvH1D*)f1.Get(Form("h_%s_mc_nobck_unfold_effcor_cross_section",varspec.c_str()));//3
  MnvH1D* mcMnv_qelike_qe   = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_qe",varspec.c_str()));//Get from N track
  MnvH1D* mcMnv_qelike_res  = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_res",varspec.c_str()));//Get from N track
  MnvH1D* mcMnv_qelike_dis  = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_dis",varspec.c_str()));//Get from N track
  MnvH1D* mcMnv_qelike_2p2h = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_2p2h",varspec.c_str()));//Get from N track
  MnvH1D* mcMnv_qelike_2p2h_nolowrecoil = (MnvH1D*)f2.Get(Form("h_%s_cross_section_qelike_2p2h",varspec.c_str()));//Get from N track
  MnvH2D* hisy_temp_mnv     = (MnvH2D*)f1.Get(Form("h_%s_qelike_templatebinning",varspec.c_str()));//Get from N track

  //Originall binning used to determine encoded to projected slices
  TH2* hisy_temp      = new TH2D(hisy_temp_mnv->GetCVHistoWithStatError());
 

  TH1* dataStat       = new TH1D(dataMnv->GetCVHistoWithStatError());
  TH1* data           = new TH1D(dataMnv->GetCVHistoWithError());
  TH1* mc             = new TH1D(mcMnv->GetCVHistoWithStatError());
  TH1* nomGenie       = new TH1D(nomGenieMnv->GetCVHistoWithStatError());
  TH1* bestGenie      = new TH1D(bestGenieMnv->GetCVHistoWithStatError());
  TH1* mc_qelike_qe   = new TH1D(mcMnv_qelike_qe->GetCVHistoWithStatError());
  TH1* mc_qelike_res  = new TH1D(mcMnv_qelike_res->GetCVHistoWithStatError());
  TH1* mc_qelike_dis  = new TH1D(mcMnv_qelike_dis->GetCVHistoWithStatError());
  TH1* mc_qelike_2p2h = new TH1D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());
  TH1* mc_qelike_2p2h_nolowrecoil = new TH1D(mcMnv_qelike_2p2h_nolowrecoil->GetCVHistoWithStatError());

  double scale = 1e39;
  if(doRatio){
    scale = 1;
    dataStat->Divide(mc);

    data->Divide(mc);
    nomGenie->Divide(mc);
    bestGenie->Divide(mc);
    mc_qelike_qe->Divide(mc);
    mc_qelike_res->Divide(mc);
    mc_qelike_dis->Divide(mc);
    mc_qelike_2p2h->Divide(mc);
    mc_qelike_2p2h_nolowrecoil->Divide(mc);
    mc->Divide(mc);
  }


  vector<int> mycolors  = MnvColors::GetColors(2);
  
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);

  nomGenie->SetLineColor(kBlue);
  nomGenie->SetLineWidth(2);
  bestGenie->SetLineColor(mycolors[10]);
  bestGenie->SetLineWidth(2);
  mc_qelike_qe->SetLineColor(mycolors[3]);
  mc_qelike_res->SetLineColor(mycolors[4]);
  mc_qelike_dis->SetLineColor(mycolors[5]);
  mc_qelike_2p2h->SetLineColor(mycolors[16]);
  mc_qelike_2p2h_nolowrecoil->SetLineColor(mycolors[16]);
  mc_qelike_2p2h_nolowrecoil->SetLineStyle(2);
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.6);
  //data->SetMarkerStyle(1);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);
  
  dataStat->SetLineColor(kBlack);
  dataStat->SetMarkerStyle(1);


  //Convert encoded space into the N 1D slices
  vector <TH1*> dataStat1D       = PrepareHist(hisy_temp, dataStat, "dataStat1D",scale,doRatio);
  vector <TH1*> data1D           = PrepareHist(hisy_temp, data, "data1D",scale,doRatio);
  vector <TH1*> mc1D             = PrepareHist(hisy_temp, mc, "mc1D",scale,doRatio);
  vector <TH1*> nomGenie1D       = PrepareHist(hisy_temp, nomGenie, "nomGenie1D",scale,doRatio);
  vector <TH1*> bestGenie1D      = PrepareHist(hisy_temp, bestGenie, "bestGenie1D",scale,doRatio);
  vector <TH1*> mc_qelike_qe1D   = PrepareHist(hisy_temp, mc_qelike_qe, "mc_qelike_qe1D",scale,doRatio);
  vector <TH1*> mc_qelike_res1D  = PrepareHist(hisy_temp, mc_qelike_res, "mc_qelike_res1D",scale,doRatio);
  vector <TH1*> mc_qelike_dis1D  = PrepareHist(hisy_temp, mc_qelike_dis, "mc_qelike_dis1D",scale,doRatio);
  vector <TH1*> mc_qelike_2p2h1D = PrepareHist(hisy_temp, mc_qelike_2p2h, "mc_qelike_2p2h1D",scale,doRatio);
  vector <TH1*> mc_qelike_2p2h_nolowrecoil1D = PrepareHist(hisy_temp, mc_qelike_2p2h_nolowrecoil, "mc_qelike_2p2h1D_nolowrecoil",scale,doRatio);
  
  //Now all the histograms classes.
  vector <std::pair<vector <TH1*>, const char*>> allHists;
  allHists.push_back(std::make_pair(dataStat1D,"histep"));
  allHists.push_back(std::make_pair(mc1D,"histl"));
  if(doGenies){  
    allHists.push_back(std::make_pair(nomGenie1D,"histl"));
    allHists.push_back(std::make_pair(bestGenie1D,"histl"));
  } else {
    allHists.push_back(std::make_pair(mc_qelike_qe1D,"histl"));
    allHists.push_back(std::make_pair(mc_qelike_res1D,"histl"));
    allHists.push_back(std::make_pair(mc_qelike_dis1D,"histl"));
    allHists.push_back(std::make_pair(mc_qelike_2p2h1D,"histl"));
    allHists.push_back(std::make_pair(mc_qelike_2p2h_nolowrecoil1D,"histl"));
  }
  allHists.push_back(std::make_pair(data1D,"histep"));

  

  vector<double> multi = GetScales(allHists,6.99,0.75);


  GridCanvas* gc=NULL;//new GridCanvas(uniq(), grid_x, grid_y, 1000, 700);
  gc = plotNonGrid1D(allHists,varLatex.at(varx),varLatex.at(vary),hisy_temp_mnv,doMultipliers ? &multi[0]:NULL);
  string ytitle = "d^{2}#sigma/"+varLatex.at(varx)+varLatex.at(vary)+" (x10^{-39} cm^{2}/GeV^{2}/c^{2}/Nucleon)";
  if(doRatio)gc->SetYTitle("Ratio to MINERvA Tune v1.0.1");
  else gc->SetYTitle(ytitle.c_str());
  gc->Modified();
  
  TLegend* leg=new TLegend(0.82, 0.1, 0.9, 0.35);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MnvGENIE", "l");
  if(doGenies){
    leg->AddEntry(nomGenie, "GENIE 2.8.4", "l");
    leg->AddEntry(bestGenie, "2p2h and RPA", "l");
  }
  else{
    leg->AddEntry(mc_qelike_qe,"QELike-QE","l");
    leg->AddEntry(mc_qelike_res,"QELike-Resonant","l");
    leg->AddEntry(mc_qelike_dis,"QELike-DIS","l");
    leg->AddEntry(mc_qelike_2p2h,"QELike-2p2h","l");
    leg->AddEntry(mc_qelike_2p2h_nolowrecoil,"2p2h without fit","l");
    
  }
  leg->Draw();
  if(doRatio)gc->SetYLimits(0, 1.99);
  else gc->SetYLimits(0, 6.99);
  gc->Modified();


  //  gc->Print("test.png");

  if(doGenies){
    if(doRatio){
      gc->Print(doMultipliers ? Form("nu-2d-xsec-genies-nongrid-iter-%s-%s-multiplier-ratio.png",iiter.c_str(), var.c_str()) : Form("nu-2d-xsec-genies-nongrid-iter-%s-%s-ratio.png",iiter.c_str(), var.c_str()));
    }
    else{
      gc->Print(doMultipliers ? Form("nu-2d-xsec-genies-nongrid-iter-%s-%s-multiplier.png",iiter.c_str(), var.c_str()) : Form("nu-2d-xsec-genies-nongrid-iter-%s-%s.png",iiter.c_str(), var.c_str()));
    }
  }
  else {
    if(doRatio){
      gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-nongrid-iter-%s-%s-multiplier-ratio.png",iiter.c_str(), var.c_str()) : Form("nu-2d-xsec-comps-nongrid-iter-%s-%s-ratio.png",iiter.c_str(), var.c_str()));
    }
    else{
      gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-nongrid-iter-%s-%s-multiplier.png",iiter.c_str(), var.c_str()) : Form("nu-2d-xsec-comps-nongrid-iter-%s-%s.png",iiter.c_str(), var.c_str()));
    }

  }


}

int main(int argc, char* argv[])
{

  InitVarLatex();
  string a = argv[4];
  makePlots(true,false,false,argv[1],argv[2],argv[3],a);
  makePlots(false,false,false,argv[1],argv[2],argv[3],a);
  makePlots(false,false,true,argv[1],argv[2],argv[3],a);
  return 0;
}
