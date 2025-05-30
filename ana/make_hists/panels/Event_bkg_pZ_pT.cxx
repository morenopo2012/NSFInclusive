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
#include "PlotUtils/MnvColors.h"
//#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"
#include <algorithm>
using namespace PlotUtils;
using namespace std;
void makePlots(bool doMultipliers,bool doRatio,string location)
{
 // ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  gStyle->SetLabelSize(0.04);


// TFile f1("/minerva/data/users/sakhter/root_files_Feb_7/rootfiles_MADP3/Event2/Hists_EventSelection_MAD_P3_2D_Combine_t1_z26_AntiNu.root");

//TFile f1("../Projected2D.root");
//TFile f1("Bkg_Sub_Tracker_fullp3_x_Q2_2D_proj.root");
//TFile f1("ES_3D_Project_2D_fullRHC_1_26.root");
TFile f1("/exp/minerva/app/users/omorenop/cmtuser/git-Mat/NSFNukeCCInclusive/ana/make_hists/Hists_EventSelection_t99_z99_AntiNu_minervame5A.root");
//  MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_data_pZmu_x_y_zy");
//  MnvH2D* mcMnv=(MnvH2D*)f1.Get("h_mc_pZmu_x_y_zy");

  //MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_data_x_Q2_W_yx");
  MnvH2D* mcMnv=(MnvH2D*)f1.Get("h_mc_reco_Lead_Eavailable_q3");
  //MnvH2D* databkg=(MnvH2D*)f1.Get("h_data_bkg_POT_norm_yx");
  MnvH2D* mcbkg=(MnvH2D*)f1.Get("h_mc_reco_Carbon_Eavailable_q3");


  double DataPOT = 1.11707110844e+21;
  double MCPOT = 4.96197407132e+21;
   double scale = DataPOT/MCPOT;
    //dataMnv->Scale(1e-4, "width");
    mcMnv->Scale(scale*1e-4, "width");
    //databkg->Scale(1e-4,"width");
     mcbkg->Scale(scale*1e-4, "width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
//  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  //TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mcTotalError=new TH2D(mcMnv->GetCVHistoWithError());
  TH2* mc_bkg=new TH2D(mcbkg->GetCVHistoWithError());
   //TH2* data_bkg=new TH2D(databkg->GetCVHistoWithError());
  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = MnvColors::GetColors(2);
if(!doRatio){
  mc->SetLineColor(kBlack);
 mcTotalError->SetFillColor(kRed-9);
 mc->SetLineWidth(1.8);
mc->SetFillColor(kCyan-6);
 mc_bkg->SetLineColor(kGray+2);
mc_bkg->SetLineWidth(1.7);
mc_bkg->SetFillStyle(3002);
mc_bkg->SetFillColor(kYellow-9);
 //data->SetMarkerStyle(kFullCircle);
 //data->SetMarkerSize(0.5);
  //data->SetLineColor(kBlack);
    //data->SetLineWidth(2);
    //data->SetMarkerColor(kBlack);
    //  data_bkg->SetLineColor(kPink+4);
//data_bkg->SetMarkerStyle(kFullCircle);
//data_bkg->SetMarkerSize(0.5);
//data_bkg->SetMarkerColor(kPink+4);
//data_bkg->SetFillColor(kRed-10);
//data_bkg->SetFillStyle(3002);
}
else {
mc->SetLineColor(kRed);
 mcTotalError->SetFillColor(kRed-9);
 mc->SetLineWidth(2);
}

//if(doRatio){
 //data->Divide(mc);
	  //mcTotalError->Divide(mc);
  // mc->Divide(mc);

//  }


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
 if(!doRatio){

	 histAndOpts.push_back(std::make_pair(mcTotalError,       "E2"));
  //  histAndOpts.push_back(std::make_pair(mcStat,       "graphe3"));
  histAndOpts.push_back(std::make_pair(mc,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_bkg,       "hist"));
 //histAndOpts.push_back(std::make_pair(data,     "histpe1"));
  //histAndOpts.push_back(std::make_pair(data_bkg,     "histpe1"));
}
else{
histAndOpts.push_back(std::make_pair(mc,       "hist"));
 //histAndOpts.push_back(std::make_pair(data,     "histpe1"));
//1 histAndOpts.push_back(std::make_pair(dataStat, "X0" "histpe1"));
}
  // Values to multiply each bin by to get them on a similar range
double  multipliersxQ21[]={1, 3, 10,50,400,15000};
//double  multipliersxQ21[]={1, 1,1,1,1,1.5,15,80,150};


  GridCanvas* gc=plotXAxis1D(histAndOpts, "Bjorken x", "Q^{2}(GeV^{2})",3,2,850,520, doMultipliers ? multipliersxQ21 :NULL);

  // Set the y range manually. Can also use gc->Remax() to guess automatically
   gc->Remax();
//  gc->SetYLimits(0, 6000);
//  if(doRatio) gc->SetYLimits(0,1.99);
//  else  gc->SetYLimits(0, 4.59);
  if(doRatio){ gc->SetYTitle("Ratio");
     gc->SetYLimits(0,2);}
// else gc->SetYTitle("Events  per (GeV)^{2}");
    else gc->SetYTitle("Event Rate x 10^{-4}");

  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.0, 0.0, 0.4, 0.1);
TLegend* leg0=new TLegend(0.68, 0.0, 0.95, 0.1);
 TLatex* latex = new TLatex(0, 0, "#color[2]{#font[72]{MINERvA Work in Progress}}");
    leg0->AddEntry(latex, "", "");
    leg0->SetBorderSize(0);
    leg0->SetTextSize(0.035);

  leg->SetNColumns(2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.025);
  leg->SetTextFont(132);
  //leg->AddEntry(data, "MINERvA Data", "lpe1");
if(!doRatio){  leg->AddEntry(mc, "Total MC (MINERvA tuneV1)", "f");}
else{leg->AddEntry(mc, "Total MC (MINERvAtuneV1)", "l");}
if(!doRatio){
  leg->AddEntry(mc_bkg,"Total MC bkg","f");
  //leg->AddEntry(data_bkg,"Data constrained bkg prediction","lpe1");
}

  TLegend* leg2 = new TLegend(0.7, 0.0, 1, 0.1);
  leg2->SetNColumns(2);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.025);


//leg2->AddEntry(mc_dis_sis,"Wrong sign","l");
//leg2->AddEntry(mc_2p2h,"Other","l");
  TLegend* leg3=new TLegend(0.11, 0.81, 0.4, 0.91);
  leg3->SetNColumns(1);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.03);
//  leg3->AddEntry((TObject*)0, "MINERvA Preliminary", "");
  leg3->SetTextColor(4);
//  leg3->AddEntry((TObject*)0,"DataPOT:1.61e+20","");
//  leg3->AddEntry((TObject*)0,"MCPOT:6.4e+20", "");
   leg3->AddEntry((TObject*)0,"DataPOT:1.118e+21","");
  leg3->AddEntry((TObject*)0,"MCPOT:4.962e+21", "");


  TLatex *mytex = new TLatex();
  TLine *myline = new TLine();
  mytex->SetTextFont(42);
  mytex->SetTextSize(0.050);
  leg->Draw("SAME");
  leg0->Draw("SAME");
  leg2->Draw("SAME");
   leg3->Draw("SAME");
  //  mytex->DrawLatex(0.67,0.23,"GENIE Components:");
  //  myline->DrawLine(0.67,0.225,0.847,0.225);

  if(doRatio){
    gc->Print(doMultipliers ? "Event_x_Q2_target_bkg_Sub.eps" : "Event_x_Q2_target_bkg_Sub_ratio.eps");
    gc->Print(doMultipliers ? "Event_x_Q2_target_multiplier_bkg_Sub_ratio.png" : "Event_x_Q2_target_bkg_Sub_ratio.png");
//    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier_ratio.C" : "nu-2d-xsec-comps-pt_ratio.C");
   // gc->SetYLimits(0,2);

 }
  else{
    gc->Print(doMultipliers ? "Event_x_Q2_target_multiplier_bkg_Sub.eps" : "Event_x_Q2_target_bkg_Sub.eps");
    gc->Print(doMultipliers ? "Event_x_Q2_target_multiplier_bkg_Sub.png" : "Event_x_Q2_target_bkg_Sub.png");
    gc->Print(doMultipliers ? "Event_x_Q2_target_multiplier_bkg_Sub.C" : "Event_x_Q2_target_bkg_Sub.C");
  }

}

int main(int argc, char* argv[])
{
  //  makePlots(true,true,argv[1]); //multipliers on ratio don't make sense
  makePlots(true,false,argv[1]);
  makePlots(false,true,argv[1]);
  makePlots(false,false,argv[1]);
  return 0;
}
