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

void makePlots(bool doMultipliers,bool doGenies,string location)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  //TFile f1("/minerva/data/users/schellma/PAPER9/Modified_pzmu_ptmu_result_Feb2020_h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section_cv_bands_binwidth_TH.root");//Final result
  //  MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section_CV_WithStatErr");

  TFile f1("/minerva/data/users/schellma/PAPER9/Modified_pzmu_ptmu_result_Feb2020.root");
  MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  dataStat->Scale(1e39,"width");

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



  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  /*
  double multipliers[]={1, 1, 1, 1,
                        2, 4, 6, 8,
                        10, 20, 30, 50};

  GridCanvas* gc=plotXAxis1D(histAndOpts, "P_{||}", "Pt", doMultipliers ? multipliers : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  gc->SetYLimits(0, 2.99);
  gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right


  if(doGenies){
    gc->Print(doMultipliers ? "nu-2d-xsec-genies-pt-multiplier.eps" : "nu-2d-xsec-genies-pt.eps");
    gc->Print(doMultipliers ? "nu-2d-xsec-genies-pt-multiplier.png" : "nu-2d-xsec-genies-pt.png");
    gc->Print(doMultipliers ? "nu-2d-xsec-genies-pt-multiplier.C" : "nu-2d-xsec-genies-pt.C");
  }
  else{
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.eps" : "nu-2d-xsec-comps-pt.eps");
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.png" : "nu-2d-xsec-comps-pt.png");
    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier.C" : "nu-2d-xsec-comps-pt.C");
  }
  */
  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  double  multipliers2[]={1, 1, 1, 1, 1,
			  1, 1, 1, 1, 1,
                          2, 4, 10, 30, 100};

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=plotYAxis1D(histAndOpts, "Pt", "P_{||}",5,3, 800,500,doMultipliers ? multipliers2 : NULL);
  gc2->SetYLimits(0, 2.7);
  gc2->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  gc2->Modified();
  if(doGenies){
    gc2->Print(doMultipliers ? "nu-2d-xsec-genies-pz-multiplier.eps" : "nu-2d-xsec-genies-pz.eps");
    gc2->Print(doMultipliers ? "nu-2d-xsec-genies-pz-multiplier.png" : "nu-2d-xsec-genies-pz.png");
    gc2->Print(doMultipliers ? "nu-2d-xsec-genies-pz-multiplier.C" : "nu-2d-xsec-genies-pz.C");
  }
  else{
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.eps" : "nu-2d-xsec-comps-pz.eps");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.png" : "nu-2d-xsec-comps-pz.png");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.C" : "nu-2d-xsec-comps-pz.C");
  }

}

int main(int argc, char* argv[])
{
  makePlots(true,true,argv[1]);
  makePlots(true,false,argv[1]);
  makePlots(false,true,argv[1]);
  makePlots(false,false,argv[1]);

  return 0;
}
