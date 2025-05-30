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

void makePlots(bool doMultipliers,bool doRatio,string location, string variable)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  gStyle->SetLabelSize(0.04);


  string xvar = "";
  string yvar = "";

  bool logx = false;
  bool logy = false;

  double ymax = 4.99;
  double ratiomax = 1.99;

  if(variable=="pzmu_ptmu"){
    xvar = "p_{||}";
    yvar = "p_{t}";
    
  }
  else if(variable=="hadQ2_ptmu"){
    xvar = "Q^{2}";
    yvar = "p_{t}";
    logx=true;
  }
  else if(variable=="nu_hadQ2"){
    xvar = "#nu";
    yvar = "Q^{2}";
    logy=true;
    logx=true;
    ymax=9.99;
    ratiomax = 3.99;
  }
  else if(variable=="W_hadQ2"){
    xvar = "W";
    yvar = "Q^{2}";
    logx=true;
    logy=true;
    ymax = 9.99;
    ratiomax = 4.99;
  }

  else if(variable=="hadQ2_x"){
    xvar = "Q^{2}";
    yvar = "x";
    logx=true;
    logy=true;
    ymax = 1329.99;
    ratiomax = 4.99;
  }
  else if(variable=="hadQ2_y"){
    xvar = "Q^{2}";
    yvar = "y";
    logx=true;
    logy=true;
    ymax = 59.99;
    ratiomax = 2.99;
  }
  else if(variable=="x_y"){
    xvar = "x";
    yvar = "y";
    logx=true;
    logy=true;
    ymax = 39.99;
    ratiomax = 1.99;
  }


  



  TFile f1(location.c_str());//Final result
  MnvH2D* dataMnv=(MnvH2D*)f1.Get(Form("h_%s_data_nobck_unfold_effcor_cross_section",variable.c_str()));
  MnvH2D* mcMnv=(MnvH2D*)f1.Get(Form("h_%s_mc_nobck_unfold_effcor_cross_section",variable.c_str()));


  dataMnv->GetXaxis()->SetLabelSize(0.04);
  dataMnv->GetYaxis()->SetLabelSize(0.04);


  MnvH2D* mcMnv_qe = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qe",variable.c_str()));//Get from N track
  MnvH2D* mcMnv_res = (MnvH2D*)f1.Get(Form("h_%s_cross_section_res",variable.c_str()));//Get from N track
  MnvH2D* mcMnv_dis = (MnvH2D*)f1.Get(Form("h_%s_cross_section_dis",variable.c_str()));//Get from N track
  MnvH2D* mcMnv_dis_dis = (MnvH2D*)f1.Get(Form("h_%s_cross_section_dis_dis",variable.c_str()));//Get from N track
  MnvH2D* mcMnv_dis_sis = (MnvH2D*)f1.Get(Form("h_%s_cross_section_dis_sis",variable.c_str()));//Get from N track
  MnvH2D* mcMnv_2p2h = (MnvH2D*)f1.Get(Form("h_%s_cross_section_2p2h",variable.c_str()));//Get from N track
  MnvH2D* mcMnv_oth = (MnvH2D*)f1.Get(Form("h_%s_cross_section_oth",variable.c_str()));//Get from N track

  dataMnv->GetXaxis()->SetTitle(xvar.c_str());

  dataMnv->Scale(1e39, "width");
  mcMnv->Scale(1e39, "width");

  mcMnv_qe->Scale(1e39,"width");
  mcMnv_res->Scale(1e39,"width");
  mcMnv_dis->Scale(1e39,"width");
  mcMnv_dis_dis->Scale(1e39,"width");
  mcMnv_dis_sis->Scale(1e39,"width");
  mcMnv_2p2h->Scale(1e39,"width");
  mcMnv_oth->Scale(1e39,"width");

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

  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = MnvColors::GetColors(9);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);

  mc_qe->SetLineColor(mycolors[2]);
  mc_res->SetLineColor(mycolors[0]);
  mc_dis_dis->SetLineColor(kViolet-3);
  mc_dis_sis->SetLineColor(mycolors[1]);
  mc_oth->SetLineColor(mycolors[4]);


  //Add 2p2h with qe to reduce number of categories
  mc_qe->Add(mc_2p2h);

  mc_qe->SetLineWidth(1.5);
  mc_res->SetLineWidth(1.5);
  mc_dis_dis->SetLineWidth(1.5);
  mc_dis_sis->SetLineWidth(1.5);
  mc_oth->SetLineWidth(1.5);

  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.5);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);

  dataStat->SetLineColor(kBlack);


  if(doRatio){
    //    TH2 *tmpden = (TH2*)mc->Clone("tmpden");
    //    tmpden->Sumw2(false);
    data->Divide(mc);
    dataStat->Divide(mc);
    mc_qe->Divide(mc);
    mc_res->Divide(mc);
    mc_dis_dis->Divide(mc);
    mc_dis_sis->Divide(mc);
    mc_oth->Divide(mc);
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
  histAndOpts.push_back(std::make_pair(dataStat, "histpe1"));
  histAndOpts.push_back(std::make_pair(mc,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(mc_qe,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(mc_res,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(mc_dis_dis,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(mc_dis_sis,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(mc_oth,       "graph0LX"));
  histAndOpts.push_back(std::make_pair(data,     "histpe1"));



  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers = GetScales(histAndOpts, true, ymax,0.75);

  GridCanvas* gc=plotXAxis1D(histAndOpts, xvar.c_str(), yvar.c_str(),doMultipliers ? &multipliers[0] : NULL);
  gc->SetLogx(logx);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  if(doRatio) gc->SetYLimits(0,ratiomax);
  else  gc->SetYLimits(0, ymax);
  if(doRatio) gc->SetYTitle("Ratio data/MINERvA Tune v1");
  else gc->SetYTitle(Form("d^{2}#sigma/d%sd%s (x10^{-39} cm^{2}/[UNIT]/Nucleon)",xvar.c_str(),yvar.c_str()));
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.6, 0.05, 0.95, 0.32);
  leg->SetNColumns(2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MINERvA Tune v1", "l");
  leg->AddEntry(mc_qe,"QE+2p2h","l");
  leg->AddEntry(mc_res,"Resonant","l");
  leg->AddEntry(mc_dis_dis,"True DIS","l");
  leg->AddEntry(mc_dis_sis,"Soft DIS","l");
  leg->AddEntry(mc_oth,"Other CC","l");

  leg->Draw("SAME");
  if(doRatio){
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier_ratio.eps" : "nu-2d-xsec-comps-pt_ratio.eps");
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier_ratio.png" : "nu-2d-xsec-comps-pt_ratio.png");
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier_ratio.C" : "nu-2d-xsec-comps-pt_ratio.C");
  }
  else{
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.eps" : "nu-2d-xsec-comps-pt.eps");
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.png" : "nu-2d-xsec-comps-pt.png");
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.C" : "nu-2d-xsec-comps-pt.C");
  }
  

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers2 = GetScales(histAndOpts, false, ymax,0.75);

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  //GridCanvas* gc2=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV/c)", "p_{||}",4,4,800,500,doMultipliers ? &multipliers2[0] : NULL);
  GridCanvas* gc2=plotYAxis1D(histAndOpts, yvar.c_str(), xvar.c_str(),doMultipliers ? &multipliers2[0] : NULL);
  gc2->SetLogx(logy);
  if(doRatio) gc2->SetYLimits(0,ratiomax);
  else  gc2->SetYLimits(0, ymax);
  if(doRatio) gc2->SetYTitle("Ratio data/MINERvA Tune v1");
  else gc2->SetYTitle(Form("d^{2}#sigma/d%sd%s (x10^{-39} cm^{2}/[UNIT]/Nucleon)",xvar.c_str(),yvar.c_str()));
  gc2->Modified();
  if(doRatio){
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.eps" : "nu-2d-xsec-comps-pz_ratio.eps");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.png" : "nu-2d-xsec-comps-pz_ratio.png");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.C" : "nu-2d-xsec-comps-pz_ratio.C");
  }
  else{
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.eps" : "nu-2d-xsec-comps-pz.eps");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.png" : "nu-2d-xsec-comps-pz.png");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.C" : "nu-2d-xsec-comps-pz.C");
  }

}

int main(int argc, char* argv[])
{
  //  makePlots(true,true,argv[1]); //multipliers on ratio don't make sense
  makePlots(true,false,argv[1],argv[2]);
  makePlots(false,true,argv[1],argv[2]);
  makePlots(false,false,argv[1],argv[2]);
  return 0;
}
