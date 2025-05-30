#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
// #include "../util/plot/plot.h"

#include "TGaxis.h"
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

//#include "plot_ME.h"
#include "plot.h"
#include <algorithm>
using namespace PlotUtils;

void makePlots(bool doMultipliers, bool doPrelimLabel = false)
{
  //ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  //gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 1-track, 2+track, N-track

  //TFile f1("Hists_EventSelection_t1_z26_Nu_v1_NuE.root");//1track
  TFile f1("/minerva/app/users/zdar/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive/ana/make_hists/CrossSection_NSF/XS_pzmu_t14_z82.root");//1track
  //TFile f1("../../make_hists/FullTarget/Hists_EventSelection_t1_z82_Nu_v1_.root");//1track
  

  MnvH2D* mcMnv=(MnvH2D*)f1.Get("h_recoCrossSection");
  MnvH2D* dataMnv = (MnvH2D*)f1.Get("h_dataCrossSection");

  


  dataMnv->Scale(1e-3, "width");
  mcMnv->Scale(1e-3, "width");
  mcMnv->Scale(0.467);
  



  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithError());
  TH2* mcStat=new TH2D(mcMnv->GetCVHistoWithStatError());

  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  //  mcStat->SetFillStyle(1001);
  mcStat->SetFillColor(kPink + 1);
  mc->SetLineWidth(2);
  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.7);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);

  dataStat->SetMarkerStyle(1);
  dataStat->SetLineColor(kBlack);
  dataStat->SetLineWidth(2);

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
  histAndOpts.push_back(std::make_pair(mcStat,       "e2"));
  histAndOpts.push_back(std::make_pair(mc,       "hist"));
  //histAndOpts.push_back(std::make_pair(mcStat,       "e3"));
  //histAndOpts.push_back(std::make_pair(dataStat, "histpe1"));
  histAndOpts.push_back(std::make_pair(data,     "histpe1"));



  // ----------------------------------------------------------------------------------
  //
  // First make Ehad in bins of Emu

  // Values to multiply each bin by to get them on a similar range
  double multipliers3r[]={1, 1, 1, 1,
                        1, 1, 1, 1,
                        1, 1, 1, 1};

  GridCanvas* gc=plotXfrac1D(histAndOpts, "Muon Energy (GeV)", "Ehad", doMultipliers ? multipliers3r : NULL);
  //GridCanvas* gc=plotXfrac1D(histAndOpts, "W (GeV)", "Q2", doMultipliers ? multipliers3r : NULL);
  //GridCanvas* gc=plotXfrac1D(histAndOpts, "x", "Q2", doMultipliers ? multipliers3r : NULL);
  //GridCanvas* gc=plotXfrac1D(histAndOpts, "x", "y", doMultipliers ? multipliers3r : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  //gc->Remax();
  gc->SetYLimits(0, 2);
  gc->SetYTitle("Data / MnvGENIE_v1");
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.17, 0.7, 0.31, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.023);
  leg->AddEntry(" ", "MINERvA Preliminary");
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MINERvA Tune", "l");
  //leg->AddEntry(" ","Lead of Target 5");
  leg->AddEntry(" ","Tracker");
  //leg->AddEntry(" ","Lead of Target 3");
  leg->Draw("SAME");
    gc->Print(doMultipliers ? "nu-2d-XS_ratios-Tracker-model-Ehad-multiplier.png" : "nu-2d-XS-ratios-Tracker-model-Ehad.png");
    //gc->Print(doMultipliers ? "nu-2d-evtrate-model-Q2-multiplier.png" : "nu-2d-evtrate-model-Q2.png");
    //gc->Print(doMultipliers ? "nu-2d-evtrate-model-Q2-multiplier.png" : "nu-2d-evtrate-model-Q2.png");
    //gc->Print(doMultipliers ? "nu-2d-evtrate-model-y-multiplier.png" : "nu-2d-evtrate-model-y.png");
  


  // ------------------------------------------------------------------------------
  //
  // Now make Emu in bins of Ehad. It's all the same

  // Values to multiply each bin by to get them on a similar range
  double  multipliers7[]={1,1, 1, 1, 1,
			  1, 1, 1, 1, 1,
			  1, 1, 1};
  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=plotYfrac1D(histAndOpts, "Hadronic Energy (GeV)", "Emu", doMultipliers ? multipliers7 : NULL);
  //GridCanvas* gc2=plotYfrac1D(histAndOpts, "Q2 (GeV^{2})", "W", doMultipliers ? multipliers7 : NULL);
  //GridCanvas* gc2=plotYfrac1D(histAndOpts, "Q2 (GeV^{2})", "x", doMultipliers ? multipliers7 : NULL);
  //GridCanvas* gc2=plotYfrac1D(histAndOpts, "y", "x", doMultipliers ? multipliers7 : NULL);
  //gc2->Remax();
  gc2->SetYLimits(0, 2);
  gc2->SetYTitle("Data / MnvGENIE_v1");
  gc2->Modified();
    

  gc2->Print(doMultipliers ? "nu-2d-XS-ratios-Tracker-model-Emu-multiplier.png" : "nu-2d-XS-ratios-Tracker-model-Emu.png");
  //gc2->Print(doMultipliers ? "nu-2d-evtrate-model-W-multiplier.png" : "nu-2d-evtrate-model-W.png");
  //gc2->Print(doMultipliers ? "nu-2d-evtrate-model-x-multiplier.png" : "nu-2d-evtrate-model-x.png");
  //gc2->Print(doMultipliers ? "nu-2d-evtrate-model-x-multiplier.png" : "nu-2d-evtrate-model-x.png");

}

int main(int argc, char* argv[])
{

  string location = argv[1];
  //multipliers
    makePlots(false);
  //makePlots(true,true,true,location);
  //makePlots(true,false,false,location);
  //standard
/*  makePlots(false,true,false,location);
  makePlots(false,true,true,location);
  makePlots(false,false,false,location);
*/
  return 0;
}
