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
//#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"
#include <algorithm>
using namespace PlotUtils;
vector<double> GetScales(std::vector<std::pair<TH2*, const char*> >histopts, bool pzProj){
  vector<double> tmpvect;
  int nbins = histopts[0].first->GetNbinsX()+2;
  if(pzProj) nbins = histopts[0].first->GetNbinsY()+2;
  for(int i=1;i<nbins-1;i++){
    double maxval = 0;
    for(uint j=0;j<histopts.size();j++){
      TH1D *tmp = pzProj? histopts[j].first->ProjectionX("tmp",i,i): histopts[j].first->ProjectionY("tmp",i,i);
      int maxbin = tmp->GetMaximumBin();
      double content = tmp->GetBinContent(maxbin);
      if(content>maxval) maxval=content;
    }
    //we want abaout 75% of the 1.49, so 1.15
    double scale = 150./maxval;
    if(scale>1){
      int tmpscale = floor(scale*10);
      scale = tmpscale/10.0;
    }
    else{
      int tmpscale = floor(scale*10);
      scale = tmpscale/10.0;
    }
    cout << scale << endl;
    tmpvect.push_back(scale);
  }
    return tmpvect;
}

void makePlots(bool doMultipliers, string location)
{
//  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 1-track, 2+track, N-track
  //CV
  TFile f1(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_CombinedPlaylists.root",location.c_str()));//1track
  TFile f2(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_CombinedPlaylists.root",location.c_str()));//2track

  //need pzmuptmu
  MnvH2D* mcMnv_1track=(MnvH2D*)f1.Get("h_pzmu_ptmu_mc");//Get from N track
  MnvH2D* dataMnv_1track = (MnvH2D*)f1.Get("h_pzmu_ptmu_data");//Get from N track
  MnvH2D* mcMnv_2track=(MnvH2D*)f2.Get("h_pzmu_ptmu_mc");//Get from N track
  MnvH2D* dataMnv_2track = (MnvH2D*)f2.Get("h_pzmu_ptmu_data");//Get from N track 
  MnvH2D* mcMnv_alltrack = mcMnv_1track->Clone("mc_alltrack");
  mcMnv_alltrack->Add(mcMnv_2track);
  MnvH2D* dataMnv_alltrack = dataMnv_1track->Clone("data_alltrack");
  dataMnv_alltrack->Add(dataMnv_2track);

  //bkgs
  MnvH2D* mcMnv_1track_single_neutral_pion = (MnvH2D*)f1.Get("h_pzmu_ptmu_single_neutral_pion");//Get from 1 track
  MnvH2D* mcMnv_2track_single_neutral_pion = (MnvH2D*)f2.Get("h_pzmu_ptmu_single_neutral_pion");//Get from 2+ track
  MnvH2D* mcMnv_alltrack_single_neutral_pion = mcMnv_1track_single_neutral_pion->Clone("snp_alltrack");
  mcMnv_alltrack_single_neutral_pion->Add(mcMnv_2track_single_neutral_pion);

  MnvH2D* mcMnv_1track_multi_pion = (MnvH2D*)f1.Get("h_pzmu_ptmu_multi_pion");//Get from 1 track
  MnvH2D* mcMnv_2track_multi_pion = (MnvH2D*)f2.Get("h_pzmu_ptmu_multi_pion");//Get from 2+ track
  MnvH2D* mcMnv_alltrack_multi_pion = mcMnv_1track_multi_pion->Clone("mp_alltrack");
  mcMnv_alltrack_multi_pion->Add(mcMnv_2track_multi_pion);

  MnvH2D* mcMnv_1track_oth = (MnvH2D*)f1.Get("h_pzmu_ptmu_oth");//Get from 1 track
  MnvH2D* mcMnv_2track_oth = (MnvH2D*)f2.Get("h_pzmu_ptmu_oth");//Get from 2+ track
  MnvH2D* mcMnv_alltrack_oth = mcMnv_1track_oth->Clone("mp_alltrack");
  mcMnv_alltrack_oth->Add(mcMnv_2track_oth);

  //Sig
  MnvH2D* mcMnv_1track_sig = (MnvH2D*)f1.Get("h_pzmu_ptmu_single_charged_pion");//Get from 1 track
  MnvH2D* mcMnv_2track_sig = (MnvH2D*)f2.Get("h_pzmu_ptmu_single_charged_pion");//Get from 2+ track
  MnvH2D* mcMnv_alltrack_sig = mcMnv_1track_sig->Clone("sig_alltrack");
  mcMnv_alltrack_sig->Add(mcMnv_2track_sig);


  //full rates
  dataMnv_1track->Scale(1e-3, "width");
  mcMnv_1track->Scale(1e-3, "width");
  dataMnv_2track->Scale(1e-3, "width");
  mcMnv_2track->Scale(1e-3, "width");
  dataMnv_alltrack->Scale(1e-3,"width");
  mcMnv_alltrack->Scale(1e-3, "width");

  //components
  mcMnv_1track_sig->Scale(1e-3, "width");
  mcMnv_2track_sig->Scale(1e-3, "width");
  mcMnv_alltrack_sig->Scale(1e-3, "width");

  mcMnv_1track_single_neutral_pion->Scale(1e-3, "width");
  mcMnv_2track_single_neutral_pion->Scale(1e-3, "width");
  mcMnv_alltrack_single_neutral_pion->Scale(1e-3, "width");

  mcMnv_1track_multi_pion->Scale(1e-3, "width");
  mcMnv_2track_multi_pion->Scale(1e-3, "width");
  mcMnv_alltrack_multi_pion->Scale(1e-3, "width");
  
  mcMnv_1track_oth->Scale(1e-3, "width");
  mcMnv_2track_oth->Scale(1e-3, "width");
  mcMnv_alltrack_oth->Scale(1e-3, "width");


  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat_1track=new TH2D(dataMnv_1track->GetCVHistoWithStatError());
  TH2* dataStat_2track=new TH2D(dataMnv_2track->GetCVHistoWithStatError());
  TH2* dataStat_alltrack=new TH2D(dataMnv_alltrack->GetCVHistoWithStatError());

  TH2* data_1track=new TH2D(dataMnv_1track->GetCVHistoWithError());
  TH2* data_2track=new TH2D(dataMnv_2track->GetCVHistoWithError());
  TH2* data_alltrack=new TH2D(dataMnv_alltrack->GetCVHistoWithError());

  TH2* mc_1track=new TH2D(mcMnv_1track->GetCVHistoWithStatError());
  TH2* mc_2track=new TH2D(mcMnv_2track->GetCVHistoWithStatError());
  TH2* mc_alltrack=new TH2D(mcMnv_alltrack->GetCVHistoWithStatError());


  TH2* mc_1track_sig=new TH2D(mcMnv_1track_sig->GetCVHistoWithStatError());
  TH2* mc_2track_sig=new TH2D(mcMnv_2track_sig->GetCVHistoWithStatError());
  TH2* mc_alltrack_sig=new TH2D(mcMnv_alltrack_sig->GetCVHistoWithStatError());

  TH2* mc_1track_single_neutral_pion=new TH2D(mcMnv_1track_single_neutral_pion->GetCVHistoWithStatError());
  TH2* mc_2track_single_neutral_pion=new TH2D(mcMnv_2track_single_neutral_pion->GetCVHistoWithStatError());
  TH2* mc_alltrack_single_neutral_pion=new TH2D(mcMnv_alltrack_single_neutral_pion->GetCVHistoWithStatError());

  TH2* mc_1track_multi_pion=new TH2D(mcMnv_1track_multi_pion->GetCVHistoWithStatError());
  TH2* mc_2track_multi_pion=new TH2D(mcMnv_2track_multi_pion->GetCVHistoWithStatError());
  TH2* mc_alltrack_multi_pion=new TH2D(mcMnv_alltrack_multi_pion->GetCVHistoWithStatError());

  TH2* mc_1track_oth=new TH2D(mcMnv_1track_oth->GetCVHistoWithStatError());
  TH2* mc_2track_oth=new TH2D(mcMnv_2track_oth->GetCVHistoWithStatError());
  TH2* mc_alltrack_oth=new TH2D(mcMnv_alltrack_oth->GetCVHistoWithStatError());


  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(2);
  mc_1track->SetLineColor(kRed);
  mc_1track->SetLineWidth(2);
  mc_2track->SetLineColor(kRed);
  mc_2track->SetLineWidth(2);
  mc_alltrack->SetLineColor(kRed);
  mc_alltrack->SetLineWidth(2);

  
  // These line and marker styles will be propagated to the 1D plots
  data_1track->SetMarkerStyle(kFullCircle);
  data_1track->SetMarkerSize(0.7);
  data_1track->SetLineColor(kBlack);
  data_1track->SetLineWidth(2);

  data_2track->SetMarkerStyle(kFullCircle);
  data_2track->SetMarkerSize(0.7);
  data_2track->SetLineColor(kBlack);
  data_2track->SetLineWidth(2);

  data_alltrack->SetMarkerStyle(kFullCircle);
  data_alltrack->SetMarkerSize(0.7);
  data_alltrack->SetLineColor(kBlack);
  data_alltrack->SetLineWidth(2);

  dataStat_1track->SetMarkerStyle(1);
  dataStat_1track->SetLineColor(kBlack);
  dataStat_1track->SetLineWidth(2);

  dataStat_2track->SetMarkerStyle(1);
  dataStat_2track->SetLineColor(kBlack);
  dataStat_2track->SetLineWidth(2);

  dataStat_alltrack->SetMarkerStyle(1);
  dataStat_alltrack->SetLineColor(kBlack);
  dataStat_alltrack->SetLineWidth(2);

  mc_1track_sig->SetLineColor(mycolors[5]);
  mc_2track_sig->SetLineColor(mycolors[5]);
  mc_alltrack_sig->SetLineColor(mycolors[5]);

  mc_1track_single_neutral_pion->SetLineColor(mycolors[6]);
  mc_2track_single_neutral_pion->SetLineColor(mycolors[6]);
  mc_alltrack_single_neutral_pion->SetLineColor(mycolors[6]);

  mc_1track_multi_pion->SetLineColor(mycolors[7]);
  mc_2track_multi_pion->SetLineColor(mycolors[7]);
  mc_alltrack_multi_pion->SetLineColor(mycolors[7]);

  mc_1track_oth->SetLineColor(mycolors[8]);
  mc_2track_oth->SetLineColor(mycolors[8]);
  mc_alltrack_oth->SetLineColor(mycolors[8]);



  // Make a list of the histograms we want to draw, along with the
  // draw options we want to use for them. You can add "graph" to the
  // draw options if you want the histogram to be converted to a graph
  // and then drawn. In that case the draw options are interpreted as
  // options to TGraphErrors::Draw().
  //
  // I don't know what happens if you put a "graph" first in the list,
  // so don't do that. Make sure the first item doesn't have "graph"
  // in its options
  std::vector<std::pair<TH2*, const char*> > histAndOpts_1track;
  std::vector<std::pair<TH2*, const char*> > histAndOpts_2track;
  std::vector<std::pair<TH2*, const char*> > histAndOpts_alltrack;

  histAndOpts_1track.push_back(std::make_pair(mc_1track,"histl"));
  histAndOpts_1track.push_back(std::make_pair(mc_1track_sig,"histl"));
  histAndOpts_1track.push_back(std::make_pair(mc_1track_single_neutral_pion,"histl"));
  histAndOpts_1track.push_back(std::make_pair(mc_1track_multi_pion,"histl"));
  histAndOpts_1track.push_back(std::make_pair(mc_1track_oth,"histl"));
  histAndOpts_1track.push_back(std::make_pair(dataStat_1track,"graph e"));
  histAndOpts_1track.push_back(std::make_pair(data_1track,"graph ep"));


  histAndOpts_2track.push_back(std::make_pair(mc_2track,"histl"));
  histAndOpts_2track.push_back(std::make_pair(mc_2track_sig,"histl"));
  histAndOpts_2track.push_back(std::make_pair(mc_2track_single_neutral_pion,"histl"));
  histAndOpts_2track.push_back(std::make_pair(mc_2track_multi_pion,"histl"));
  histAndOpts_2track.push_back(std::make_pair(mc_2track_oth,"histl"));
  histAndOpts_2track.push_back(std::make_pair(dataStat_2track,"graph e"));
  histAndOpts_2track.push_back(std::make_pair(data_2track,"graph ep"));


  histAndOpts_alltrack.push_back(std::make_pair(mc_alltrack,"histl"));
  histAndOpts_alltrack.push_back(std::make_pair(mc_alltrack_sig,"histl"));
  histAndOpts_alltrack.push_back(std::make_pair(mc_alltrack_single_neutral_pion,"histl"));
  histAndOpts_alltrack.push_back(std::make_pair(mc_alltrack_multi_pion,"histl"));
  histAndOpts_alltrack.push_back(std::make_pair(mc_alltrack_oth,"histl"));
  histAndOpts_alltrack.push_back(std::make_pair(dataStat_alltrack,"graph e"));
  histAndOpts_alltrack.push_back(std::make_pair(data_alltrack,"graph ep"));
  //-------------------------------------------
  //

  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.17, 0.7, 0.31, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data_1track, "MINERvA data", "lpe");
  leg->AddEntry(mc_1track, "MINERvA Tune", "l");
  leg->AddEntry(mc_1track_sig, "Single #pi^{#pm}", "l");
  leg->AddEntry(mc_1track_single_neutral_pion, "Single #pi^{0}", "l");
  leg->AddEntry(mc_1track_multi_pion, "Multi-Pion", "l");
  leg->AddEntry(mc_1track_oth, "Non-pion", "l");

  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers = GetScales( histAndOpts_alltrack, false);
  GridCanvas* gc=plotYAxis1D(histAndOpts_alltrack, "Pt", "P||", doMultipliers ? &multipliers[0] : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  gc->SetYLimits(0, 200);
  gc->SetYTitle("Event Rate(x10^{-3}) per GeV^{2}");
  gc->Modified();
 
  leg->Draw("SAME");

  gc->Print(doMultipliers ? "nu-2d-evtrate-model-alltrack--pt-multiplier.eps" : "nu-2d-evtrate-model-alltrack--pt.eps");
  gc->Print(doMultipliers ? "nu-2d-evtrate-model-alltrack--pt-multiplier.png" : "nu-2d-evtrate-model-alltrack--pt.png");
  gc->Print(doMultipliers ? "nu-2d-evtrate-model-alltrack--pt-multiplier.C" : "nu-2d-evtrate-model-alltrack--pt.C");
  
  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same
  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers2 = GetScales( histAndOpts_alltrack, true);
  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=plotXAxis1D(histAndOpts_alltrack,"Muon P||","Pt", doMultipliers ? &multipliers2[0] : NULL);
  gc2->SetYLimits(0, 200);
  gc2->SetYTitle("Event Rate(x10^{-3}) per GeV^{2}");
  gc2->Modified();
  gc2->Print(doMultipliers ? "nu-2d-evtrate-model-alltrack--pz-multiplier.eps" : "nu-2d-evtrate-model-alltrack--pz.eps");
  gc2->Print(doMultipliers ? "nu-2d-evtrate-model-alltrack--pz-multiplier.png" : "nu-2d-evtrate-model-alltrack--pz.png");
  gc2->Print(doMultipliers ? "nu-2d-evtrate-model-alltrack--pz-multiplier.C" : "nu-2d-evtrate-model-alltrack--pz.C");
  

  vector<double> multipliers3 = GetScales( histAndOpts_1track, false);
  GridCanvas* gc3=plotYAxis1D(histAndOpts_1track, "Pt", "P||", doMultipliers ? &multipliers3[0] : NULL);
  // Set the y range manually. Can also use gc3->Remax() to guess automatically
  gc3->SetYLimits(0, 200);
  gc3->SetYTitle("Event Rate(x10^{-3}) per GeV^{2}");
  gc3->Modified();
 
  leg->Draw("SAME");

  gc3->Print(doMultipliers ? "nu-2d-evtrate-model-1track--pt-multiplier.eps" : "nu-2d-evtrate-model-1track--pt.eps");
  gc3->Print(doMultipliers ? "nu-2d-evtrate-model-1track--pt-multiplier.png" : "nu-2d-evtrate-model-1track--pt.png");
  gc3->Print(doMultipliers ? "nu-2d-evtrate-model-1track--pt-multiplier.C" : "nu-2d-evtrate-model-1track--pt.C");
  
  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same
  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers4 = GetScales( histAndOpts_1track, true);
  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc4=plotXAxis1D(histAndOpts_1track,"Muon P||","Pt", doMultipliers ? &multipliers4[0] : NULL);
  gc4->SetYLimits(0, 200);
  gc4->SetYTitle("Event Rate(x10^{-3}) per GeV^{2}");
  gc4->Modified();
  gc4->Print(doMultipliers ? "nu-2d-evtrate-model-1track--pz-multiplier.eps" : "nu-2d-evtrate-model-1track--pz.eps");
  gc4->Print(doMultipliers ? "nu-2d-evtrate-model-1track--pz-multiplier.png" : "nu-2d-evtrate-model-1track--pz.png");
  gc4->Print(doMultipliers ? "nu-2d-evtrate-model-1track--pz-multiplier.C" : "nu-2d-evtrate-model-1track--pz.C");
  
  vector<double> multipliers5 = GetScales( histAndOpts_2track, false);
  GridCanvas* gc5=plotYAxis1D(histAndOpts_2track, "Pt", "P||", doMultipliers ? &multipliers3[0] : NULL);
  // Set the y range manually. Can also use gc5->Remax() to guess automatically
  gc5->SetYLimits(0, 200);
  gc5->SetYTitle("Event Rate(x10^{-3}) per GeV^{2}");
  gc5->Modified();
 
  leg->Draw("SAME");

  gc5->Print(doMultipliers ? "nu-2d-evtrate-model-2track--pt-multiplier.eps" : "nu-2d-evtrate-model-2track--pt.eps");
  gc5->Print(doMultipliers ? "nu-2d-evtrate-model-2track--pt-multiplier.png" : "nu-2d-evtrate-model-2track--pt.png");
  gc5->Print(doMultipliers ? "nu-2d-evtrate-model-2track--pt-multiplier.C" : "nu-2d-evtrate-model-2track--pt.C");
  
  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same
  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers6 = GetScales( histAndOpts_2track, true);
  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc6=plotXAxis1D(histAndOpts_2track,"Muon P||","Pt", doMultipliers ? &multipliers4[0] : NULL);
  gc6->SetYLimits(0, 200);
  gc6->SetYTitle("Event Rate(x10^{-3}) per GeV^{2}");
  gc6->Modified();
  gc6->Print(doMultipliers ? "nu-2d-evtrate-model-2track--pz-multiplier.eps" : "nu-2d-evtrate-model-2track--pz.eps");
  gc6->Print(doMultipliers ? "nu-2d-evtrate-model-2track--pz-multiplier.png" : "nu-2d-evtrate-model-2track--pz.png");
  gc6->Print(doMultipliers ? "nu-2d-evtrate-model-2track--pz-multiplier.C" : "nu-2d-evtrate-model-2track--pz.C");
  

}

int main(int argc, char* argv[])
{

  string location = argv[1];
  //multipliers
  makePlots(true,location);
  //standard
  makePlots(false,location);


  return 0;
}
