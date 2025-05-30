//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
// #include "../util/plot/plot.h"

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
#include "localColor.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;

void makePlots(bool doMultipliers, string location, string varstring,string xvar, string yvar)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 1-track, 2+track, N-track
  //CV
  //  TFile f1(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_CombinedPlaylists.root",location.c_str()));//1track
  //  TFile f2(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_CombinedPlaylists.root",location.c_str()));//2track
  TFile f3(Form("%s_CV/EffPurity_MakeFlux-1_CombinedPlaylists.root",location.c_str()));//Ntrack
  // //CV without low recoil tune
  // TFile f4(Form("%s_pion_rpa/MuonEventSelection_MakeFlux-1_Multiplicity-0_Sample-Signal_CombinedPlaylists.root",location.c_str()));//NTrack
  // //Constraint File
  // TFile f5(Form("%s_CV/SideBandFit_CombinedPlaylists.root",location.c_str()));

  //need pzmuptmu
  MnvH2D* mcMnv=(MnvH2D*)f3.Get(Form("h_%s_qelike",varstring.c_str()));//Get from N tra
  MnvH2D* mcMnv_qelike_qe = (MnvH2D*)f3.Get(Form("h_%s_qelike_qe",varstring.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_res = (MnvH2D*)f3.Get(Form("h_%s_qelike_res",varstring.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_dis = (MnvH2D*)f3.Get(Form("h_%s_qelike_dis",varstring.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_2p2h = (MnvH2D*)f3.Get(Form("h_%s_qelike_2p2h",varstring.c_str()));//Get from N track

  MnvH2D* mcMnv_Truth=(MnvH2D*)f3.Get(Form("h_%s_truth_qelike",varstring.c_str()));//Get from N tra
  MnvH2D* mcMnv_qelike_qe_Truth = (MnvH2D*)f3.Get(Form("h_%s_truth_qelike_qe",varstring.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_res_Truth = (MnvH2D*)f3.Get(Form("h_%s_truth_qelike_res",varstring.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_dis_Truth = (MnvH2D*)f3.Get(Form("h_%s_truth_qelike_dis",varstring.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_2p2h_Truth = (MnvH2D*)f3.Get(Form("h_%s_truth_qelike_2p2h",varstring.c_str()));//Get from N track
  
  mcMnv->Scale(1e-5, "width");
  mcMnv_qelike_qe->Scale(1e-5, "width");
  mcMnv_qelike_res->Scale(1e-5, "width");
  mcMnv_qelike_dis->Scale(1e-5, "width");
  mcMnv_qelike_2p2h->Scale(1e-5, "width");

  mcMnv_Truth->Scale(1e-5, "width");
  mcMnv_qelike_qe_Truth->Scale(1e-5, "width");
  mcMnv_qelike_res_Truth->Scale(1e-5, "width");
  mcMnv_qelike_dis_Truth->Scale(1e-5, "width");
  mcMnv_qelike_2p2h_Truth->Scale(1e-5, "width");
  //  mcMnv_qelike_2p2h_no_lowrec->Scale(1e-5, "width");
  //  mcMnv_bkg->Scale(1e-5, "width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mc_qelike_qe = new TH2D(mcMnv_qelike_qe->GetCVHistoWithStatError());
  TH2* mc_qelike_res = new TH2D(mcMnv_qelike_res->GetCVHistoWithStatError());
  TH2* mc_qelike_dis = new TH2D(mcMnv_qelike_dis->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h = new TH2D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());

  TH2* mc_Truth=new TH2D(mcMnv_Truth->GetCVHistoWithStatError());
  TH2* mc_qelike_qe_Truth = new TH2D(mcMnv_qelike_qe_Truth->GetCVHistoWithStatError());
  TH2* mc_qelike_res_Truth = new TH2D(mcMnv_qelike_res_Truth->GetCVHistoWithStatError());
  TH2* mc_qelike_dis_Truth = new TH2D(mcMnv_qelike_dis_Truth->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h_Truth = new TH2D(mcMnv_qelike_2p2h_Truth->GetCVHistoWithStatError());





  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);
  mc_Truth->SetLineColor(kRed);
  mc_Truth->SetLineStyle(2);

  //need to add signal and bkg colors
  mc_qelike_qe->SetLineColor(mycolors[3]);
  mc_qelike_res->SetLineColor(mycolors[4]);
  mc_qelike_dis->SetLineColor(mycolors[5]);
  mc_qelike_2p2h->SetLineColor(mycolors[16]);

  mc_qelike_qe_Truth->SetLineColor(mycolors[3]);
  mc_qelike_res_Truth->SetLineColor(mycolors[4]);
  mc_qelike_dis_Truth->SetLineColor(mycolors[5]);
  mc_qelike_2p2h_Truth->SetLineColor(mycolors[16]);

  mc_qelike_qe_Truth->SetLineStyle(2);
  mc_qelike_res_Truth->SetLineStyle(2);
  mc_qelike_dis_Truth->SetLineStyle(2);
  mc_qelike_2p2h_Truth->SetLineStyle(2);
  

  // Make a list of the histograms we want to draw, along with the
  // draw options we want to use for them. You can add "graph" to the
  // draw options if you want the histogram to be converted to a graph
  // and then drawn. In that case the draw options are interpreted as
  // options to TGraphErrors::Draw().
  //
  // I don't know what happens if you put a "graph" first in the list,
  // so don't do that. Make sure the first item doesn't have "graph"
  // in its options
  std::vector<std::pair<TH2*, const char*> > histAndOpts;

  histAndOpts.push_back(std::make_pair(mc_Truth,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_qelike_qe_Truth,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_qelike_res_Truth,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_qelike_dis_Truth,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_Truth,       "hist"));
 
  histAndOpts.push_back(std::make_pair(mc,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_qelike_qe,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_qelike_res,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_qelike_dis,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_qelike_2p2h,       "hist"));




  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  vector<double> scale = GetScales(histAndOpts, false, 75,0.75);
  GridCanvas* gc=plotYAxis1D(histAndOpts, yvar,xvar,doMultipliers ? &scale[0] : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  gc->SetYLimits(0, 75);
  gc->SetYTitle("Event Rate(x10^{-5}) per GeV^{2}");
  
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.17, 0.7, 0.31, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(mc, "MINERvA Tune", "l");
  leg->AddEntry(mc_qelike_qe,"QE","l");
  leg->AddEntry(mc_qelike_res,"Resonant","l");
  leg->AddEntry(mc_qelike_dis,"DIS","l");
  leg->AddEntry(mc_qelike_2p2h,"2p2h","l");

  //  leg->Draw("SAME");
  gc->Print(doMultipliers ? "nu-2d-effrate-model-pt-multiplier.eps" : "nu-2d-effrate-model-pt.eps");
  gc->Print(doMultipliers ? "nu-2d-effrate-model-pt-multiplier.png" : "nu-2d-effrate-model-pt.png");
  gc->Print(doMultipliers ? "nu-2d-effrate-model-pt-multiplier.C" : "nu-2d-effrate-model-pt.C");
  
  

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  vector<double> scale2 = GetScales(histAndOpts, true, 75,0.75);
  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=plotXAxis1D(histAndOpts, xvar , yvar,doMultipliers ? &scale2[0] : NULL);

  gc2->SetYLimits(0, 75);
  gc2->SetYTitle("Event Rate(x10^{-5}) per GeV^{2}");
  

  gc2->Modified();

  gc2->Print(doMultipliers ? "nu-2d-effrate-model-pz-multiplier.eps" : "nu-2d-effrate-model-pz.eps");
  gc2->Print(doMultipliers ? "nu-2d-effrate-model-pz-multiplier.png" : "nu-2d-effrate-model-pz.png");
  gc2->Print(doMultipliers ? "nu-2d-effrate-model-pz-multiplier.C" : "nu-2d-effrate-model-pz.C");


}

int main(int argc, char* argv[])
{

  string location = argv[1];
  string varstring = argv[2];
  string xvar = argv[3];
  string yvar = argv[4];
  //multipliers
  makePlots(true,location,varstring,xvar,yvar);
  //standard
  makePlots(false,location,varstring,xvar,yvar);

  return 0;
}
