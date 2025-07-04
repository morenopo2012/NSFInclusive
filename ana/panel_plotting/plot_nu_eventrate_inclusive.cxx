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
//#include "localColor.h"
#include "PlotUtils/MnvColors.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;

void makePlots(bool doMultipliers,bool doRatio, bool doBkgOnly, bool zoom, string location)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 1-track, 2+track, N-track
  //CV
  TFile f3(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-0_Sample-Signal_CombinedPlaylists.root",location.c_str()));//Ntrack


  MnvH2D* mcMnv=(MnvH2D*)f3.Get("h_pzmu_ptmu_mc");//Get from N track
  MnvH2D* dataMnv = (MnvH2D*)f3.Get("h_pzmu_ptmu_data");//Get from N track

  MnvH2D* mcMnv_qe = (MnvH2D*)f3.Get("h_pzmu_ptmu_qe");//Get from N track
  MnvH2D* mcMnv_res = (MnvH2D*)f3.Get("h_pzmu_ptmu_res");//Get from N track
  MnvH2D* mcMnv_dis = (MnvH2D*)f3.Get("h_pzmu_ptmu_dis");//Get from N track
  MnvH2D* mcMnv_dis_dis = (MnvH2D*)f3.Get("h_pzmu_ptmu_true_dis");//Get from N track
  MnvH2D* mcMnv_dis_sis = (MnvH2D*)f3.Get("h_pzmu_ptmu_sis");//Get from N track
  MnvH2D* mcMnv_2p2h = (MnvH2D*)f3.Get("h_pzmu_ptmu_2p2h");//Get from N track
  MnvH2D* mcMnv_oth = (MnvH2D*)f3.Get("h_pzmu_ptmu_oth");//Get from N track
  MnvH2D* mcMnv_NC = (MnvH2D*)f3.Get("h_pzmu_ptmu_NC");//Get from N track
  MnvH2D* mcMnv_NC_Bkg = (MnvH2D*)f3.Get("h_pzmu_ptmu_NC_Bkg");//Get from N track
  MnvH2D* mcMnv_NC_WrongSign = (MnvH2D*)f3.Get("h_pzmu_ptmu_Bkg_Wrong_Sign");//Get from N track
  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV/c)");
  dataMnv->Scale(1e-5, "width");

  mcMnv->Scale(1e-5, "width");
  mcMnv_qe->Scale(1e-5, "width");
  mcMnv_res->Scale(1e-5, "width");
  mcMnv_dis->Scale(1e-5, "width");
  mcMnv_dis_dis->Scale(1e-5, "width");
  mcMnv_dis_sis->Scale(1e-5, "width");
  mcMnv_2p2h->Scale(1e-5, "width");
  mcMnv_oth->Scale(1e-5, "width");
  mcMnv_NC->Scale(1e-5, "width");
  mcMnv_NC_Bkg->Scale(1e-5, "width");
  mcMnv_NC_WrongSign->Scale(1e-5, "width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mc_qe = new TH2D(mcMnv_qe->GetCVHistoWithStatError());
  TH2* mc_res = new TH2D(mcMnv_res->GetCVHistoWithStatError());
  TH2* mc_dis = new TH2D(mcMnv_dis->GetCVHistoWithStatError());
  TH2* mc_dis_dis = new TH2D(mcMnv_dis_dis->GetCVHistoWithStatError());
  TH2* mc_dis_sis = new TH2D(mcMnv_dis_sis->GetCVHistoWithStatError());
  TH2* mc_2p2h = new TH2D(mcMnv_2p2h->GetCVHistoWithStatError());
  TH2* mc_oth = new TH2D(mcMnv_oth->GetCVHistoWithStatError());
  TH2* mc_NC = new TH2D(mcMnv_NC->GetCVHistoWithStatError());
  TH2* mc_NC_Bkg = new TH2D(mcMnv_NC_Bkg->GetCVHistoWithStatError());
  TH2* mc_NC_WrongSign = new TH2D(mcMnv_NC_WrongSign->GetCVHistoWithStatError());

  //Add 2p2h and QE together
  mc_qe->Add(mc_2p2h);


  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = MnvColors::GetColors(9);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);
  //need to add signal and bkg colors
  mc_qe->SetLineColor(mycolors[2]);
  mc_res->SetLineColor(mycolors[0]);
  mc_dis_dis->SetLineColor(kViolet-3);
  mc_dis_sis->SetLineColor(mycolors[1]);
  mc_oth->SetLineColor(mycolors[4]);
  mc_NC->SetLineColor(kBlack);
  mc_NC_Bkg->SetLineColor(kBlack);
  mc_NC_WrongSign->SetLineColor(kBlack);


  mc_qe->SetLineWidth(1.5);
  mc_res->SetLineWidth(1.5);
  mc_dis_dis->SetLineWidth(1.5);
  mc_dis_sis->SetLineWidth(1.5);
  mc_oth->SetLineWidth(1.5);
  mc_NC->SetLineWidth(1.5);
  mc_NC_Bkg->SetLineWidth(1.5);
  mc_NC_WrongSign->SetLineWidth(1.5);


  mc_NC_Bkg->SetLineStyle(7);
  mc_NC_WrongSign->SetLineStyle(3);
  
  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.7);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);

  dataStat->SetMarkerStyle(1);
  dataStat->SetLineColor(kBlack);
  dataStat->SetLineWidth(2);

  if(doRatio){
    data->Divide(mc);
    dataStat->Divide(mc);
    mc_qe->Divide(mc);
    mc_res->Divide(mc);
    mc_dis_dis->Divide(mc);
    mc_dis_sis->Divide(mc);
    mc_oth->Divide(mc);
    mc_NC->Divide(mc);
    mc_NC_Bkg->Divide(mc);
    mc_NC_WrongSign->Divide(mc);
    mc->Divide(mc);
  }

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

  if(!doBkgOnly){
    histAndOpts.push_back(std::make_pair(dataStat, "histp"));  
    histAndOpts.push_back(std::make_pair(mc,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(mc_qe,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(mc_res,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(mc_dis_dis,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(mc_dis_sis,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(mc_oth,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(mc_NC,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(data,     "graph ep"));
  }
  histAndOpts.push_back(std::make_pair(mc_NC,       "histl"));
  if(doBkgOnly){
    histAndOpts.push_back(std::make_pair(mc_NC_Bkg,       "histl"));
    histAndOpts.push_back(std::make_pair(mc_NC_WrongSign,       "histl"));
  }
  





  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range

  vector<double> multipliers = GetScales(histAndOpts, true, 9.99,0.75);

  GridCanvas* gc=plotXAxis1DRebinPz(histAndOpts, "Muon Longitudinal Momentum (GeV/c)", "p_{t}", 4,4,800,500,doMultipliers ? &multipliers[0] : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  if(doRatio){
    gc->SetYLimits(0, 1.99);
    gc->SetYTitle("Ratio data/MINERvA Tune v1");
    if(doBkgOnly && !zoom){
      gc->SetYLimits(0,0.099);
      gc->SetYTitle("Background Fraction");
    }
    else if(doBkgOnly && zoom){
      gc->SetYLimits(0,0.0199);
      gc->SetYTitle("Background Fraction");
    }
  }
  else{
    if(doMultipliers)gc->SetYLimits(0, 9.99);
    else gc->SetYLimits(0, 14.99);
    gc->SetYTitle("Events (x10^{-5}) per (GeV/c)^{2}");
  }
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.60, 0.05, 0.95, 0.32);
  leg->SetNColumns(2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  if(!doBkgOnly){
    leg->AddEntry(mc, "MINERvA Tune v1", "l");
    leg->AddEntry(mc_qe,"QE+2p2h","l");
    leg->AddEntry(mc_res,"Resonant","l");
    leg->AddEntry(mc_dis_dis,"True DIS","l");
    leg->AddEntry(mc_dis_sis,"Soft DIS","l");
    leg->AddEntry(mc_oth,"Other CC","l");
  }
  leg->AddEntry(mc_NC,"Background","l");
  if(doBkgOnly){
    leg->AddEntry(mc_NC_Bkg,"Background: NC","l");
    leg->AddEntry(mc_NC_WrongSign,"Background: WrongSign","l");
  }
  leg->Draw("SAME");

  TLatex *prelim = new TLatex(0.11,0.955,"MINERvA Preliminary   POT = 10.61 #times 10^{20}");
  prelim->SetNDC();
  prelim->SetTextFont(42);
  prelim->SetTextSize(0.03);
  //  prelim->Draw("SAME");
  

  if(doRatio){
    if(doBkgOnly && !zoom){
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_ratio_bkgonly.eps" : "nu-2d-evtrate-model-pt_ratio_bkgonly.eps");
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_ratio_bkgonly.png" : "nu-2d-evtrate-model-pt_ratio_bkgonly.png");
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_ratio_bkgonly.C" : "nu-2d-evtrate-model-pt_ratio_bkgonly.C");
    }
    else if(doBkgOnly && zoom){
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_ratio_bkgonly_zoomed.eps" : "nu-2d-evtrate-model-pt_ratio_bkgonly_zoomed.eps");
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_ratio_bkgonly_zoomed.png" : "nu-2d-evtrate-model-pt_ratio_bkgonly_zoomed.png");
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_ratio_bkgonly_zoomed.C" : "nu-2d-evtrate-model-pt_ratio_bkgonly_zoomed.C");
    }
    else{
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_ratio.eps" : "nu-2d-evtrate-model-pt_ratio.eps");
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_ratio.png" : "nu-2d-evtrate-model-pt_ratio.png");
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_ratio.C" : "nu-2d-evtrate-model-pt_ratio.C");
    }
  }
  else{
    if(doBkgOnly){
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_bkgonly.eps" : "nu-2d-evtrate-model-pt_bkgonly.eps");
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_bkgonly.png" : "nu-2d-evtrate-model-pt_bkgonly.png");
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier_bkgonly.C" : "nu-2d-evtrate-model-pt_bkgonly.C");
    }
    else{
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier.eps" : "nu-2d-evtrate-model-pt.eps");
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier.png" : "nu-2d-evtrate-model-pt.png");
      gc->Print(doMultipliers ? "nu-2d-evtrate-model-pt-multiplier.C" : "nu-2d-evtrate-model-pt.C");
    }
  }
  

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers2 = GetScales(histAndOpts, false, 9.99,0.75);

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  //GridCanvas* gc2=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV/c)", "p_{||}",4,4,800,500,doMultipliers ? &multipliers2[0] : NULL);
  GridCanvas* gc2=plotYAxis1DRebinPt(histAndOpts, "Muon Transverse Momentum (GeV/c)", "p_{||}",4,4,800,500,doMultipliers ? &multipliers2[0] : NULL);
  if(doRatio){
    gc2->SetYLimits(0, 1.99);
    gc2->SetYTitle("Ratio data/MINERvA Tune v1");
    if(doBkgOnly && !zoom){
      gc2->SetYLimits(0,0.099);
      gc2->SetYTitle("Background Fraction");
    }
    else if(doBkgOnly && zoom){
      gc2->SetYLimits(0,0.0199);
      gc2->SetYTitle("Background Fraction");
    }

  }
  else{
    if(doMultipliers)gc2->SetYLimits(0, 9.99);
    else gc2->SetYLimits(0, 14.99);
    gc2->SetYTitle("Events (x10^{-5}) per (GeV/c)^{2}");    
  }
  gc2->Modified();
  //  prelim->Draw("SAME");
  if(doRatio){
    if(doBkgOnly && !zoom){
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_ratio_bkgonly.eps" : "nu-2d-evtrate-model-pz_ratio_bkgonly.eps");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_ratio_bkgonly.png" : "nu-2d-evtrate-model-pz_ratio_bkgonly.png");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_ratio_bkgonly.C" : "nu-2d-evtrate-model-pz_ratio_bkgonly.C");
    }
    else if(doBkgOnly && zoom){
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_ratio_bkgonly_zoomed.eps" : "nu-2d-evtrate-model-pz_ratio_bkgonly_zoomed.eps");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_ratio_bkgonly_zoomed.png" : "nu-2d-evtrate-model-pz_ratio_bkgonly_zoomed.png");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_ratio_bkgonly_zoomed.C" : "nu-2d-evtrate-model-pz_ratio_bkgonly_zoomed.C");
    }
    else{
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_ratio.eps" : "nu-2d-evtrate-model-pz_ratio.eps");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_ratio.png" : "nu-2d-evtrate-model-pz_ratio.png");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_ratio.C" : "nu-2d-evtrate-model-pz_ratio.C");
    }
  }
  else{
    if(doBkgOnly){
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_bkgonly.eps" : "nu-2d-evtrate-model-pz_bkgonly.eps");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_bkgonly.png" : "nu-2d-evtrate-model-pz_bkgonly.png");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier_bkgonly.C" : "nu-2d-evtrate-model-pz_bkgonly.C");
    }
    else{
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier.eps" : "nu-2d-evtrate-model-pz.eps");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier.png" : "nu-2d-evtrate-model-pz.png");
      gc2->Print(doMultipliers ? "nu-2d-evtrate-model-pz-multiplier.C" : "nu-2d-evtrate-model-pz.C");
    }
  }
}

int main(int argc, char* argv[])
{

  string location = argv[1];

  makePlots(false,true,false,false,location);
  makePlots(true,false,false,false,location);
  makePlots(true,false,true,false,location);
  makePlots(false,true,true,false,location);
  makePlots(false,true,true,true,location);
  makePlots(false,false,false,false,location);

  return 0;
}
