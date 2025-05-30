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

//#include "plot_ME.h"
#include "plot.h"
#include <algorithm>
using namespace PlotUtils;

void makePlots(bool doMultipliers)
{
  //ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  //gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  //three files 1-track, 2+track, N-track

  //TFile f1("Hists_EventSelection_t1_z26_Nu_v1_NuE.root");//1track
  //TFile f1("histsEventsDIS/Hists_EventSelection_t5_z82_Nu_v1_.root");//1track
  //TFile f1("histsFullStatDIS/Hists_EventSelection_t1_z26_Nu_v1_.root");//1track
  ///////////TFile f1("histsFullStatDIS/Hists_EventSelection_t3_z06_Nu_v1_.root");//1track
  //TFile f1("/minerva/data/users/zdar/MADHists/histsFullStatDIS/histsFullStatDIS/Hists_EventSelection_t2_z26_Nu_v1_.root");//1track
  TFile f1("/minerva/app/users/zdar/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive/ana/make_hists/Background_macros/Hists_PhysicsSideband_t1_z26_Nu.root");//1track
  //need pzmuptmu
/*  MnvH2D* mcMnv=(MnvH2D*)f1.Get("h_mc_Emu_Ehad");//Get from N track
  MnvH2D* mcMnvUS=(MnvH2D*)f1.Get("h_mc_USCH_Emu_Ehad");//Get from N track
  MnvH2D* mcMnvDS=(MnvH2D*)f1.Get("h_mc_DSCH_Emu_Ehad");//Get from N track
  MnvH2D* mcMnvTrans=(MnvH2D*)f1.Get("h_mc_trans_SB_Emu_Ehad");//Get from N track
  MnvH2D* mcMnvContin=(MnvH2D*)f1.Get("h_mc_contin_SB_Emu_Ehad");//Get from N track
  MnvH2D* dataMnv = (MnvH2D*)f1.Get("h_data_Emu_Ehad");//Get from N track
 */
  MnvH2D* TinT=(MnvH2D*)f1.Get("selected_mc_trans_in_trans_Emu_Ehad");//Get from N track
  MnvH2D* CinT=(MnvH2D*)f1.Get("selected_mc_contin_in_trans_Emu_Ehad");//Get from N track
  MnvH2D* SinT=(MnvH2D*)f1.Get("selected_mc_signal_in_trans_Emu_Ehad");//Get from N track

  MnvH2D* TUS=(MnvH2D*)f1.Get("selected_trans_CH_US_Emu_Ehad");//Get from N track
  MnvH2D* TDS=(MnvH2D*)f1.Get("selected_trans_CH_DS_Emu_Ehad");//Get from N track

  
  MnvH2D* TData=(MnvH2D*)f1.Get("selected_data_trans_Emu_Ehad");//Get from N track







  //dataMnv->GetXaxis()->SetTitle("Muon Energy (GeV)");
  //dataMnv->SetTitle("Muon Energy (GeV)");

  TData->Scale(1e-3, "width");

  TinT->Scale(1e-3, "width");
  CinT->Scale(1e-3, "width");
  SinT->Scale(1e-3, "width");

  TUS->Scale(1e-3, "width");
  TDS->Scale(1e-3, "width");



  TinT->Scale(0.421);
  CinT->Scale(0.421);
  SinT->Scale(0.421);

  TUS->Scale(0.421);
  TDS->Scale(0.421);

  



  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(TData->GetCVHistoWithStatError());
  TH2* data=new TH2D(TData->GetCVHistoWithError());
  TH2* TinTrans=new TH2D(TinT->GetCVHistoWithError());
  TH2* CinTrans=new TH2D(CinT->GetCVHistoWithError());
  TH2* SinTrans=new TH2D(SinT->GetCVHistoWithError());
  TH2* TransUS=new TH2D(TUS->GetCVHistoWithError());
  TH2* TransDS=new TH2D(TDS->GetCVHistoWithError());

  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(2);
  TinTrans->SetLineColor(kRed);
  //mcStat->SetFillColorAlpha(kPink + 1, 0.4);
 // ->SetFillColor(kRed-9);
  TinTrans->SetLineWidth(2);
  CinTrans->SetLineWidth(2);
  CinTrans->SetLineColor(kBlue);
  //mcDS->SetLineColor(kGreen);
  SinTrans->SetLineWidth(2);
  SinTrans->SetLineColor(kRed);
  //mcContin->SetLineColor(kOrange);
  //mcDS->SetLineWidth(2);
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
  //histAndOpts.push_back(std::make_pair(mcStat,       "E2"));
  histAndOpts.push_back(std::make_pair(TinTrans,       "hists"));
  histAndOpts.push_back(std::make_pair(CinTrans,       "hists"));
  histAndOpts.push_back(std::make_pair(SinTrans,       "hists"));
  histAndOpts.push_back(std::make_pair(TransUS,       "hists"));
  histAndOpts.push_back(std::make_pair(TransDS,       "hists"));
  histAndOpts.push_back(std::make_pair(dataStat, "graph e"));
  histAndOpts.push_back(std::make_pair(data,     "graph ep"));



  // ----------------------------------------------------------------------------------
  //
  // First make Ehad in bins of Emu

  // Values to multiply each bin by to get them on a similar range
  //double multipliers1[]={0.4, 0.3, 0.7, 1,
  //                      3, 4,8, 2,
  //                      2, 1, 1, 1};
  double multipliers1[]={1, 60, 400, 1000,
                        1500, 2000, 3000, 300,
                        2, 1, 1, 1};
  double multipliersxQ21[]={1, 3, 8, 25,
                        100, 1, 1, 1,
                        1, 1, 1, 1};
  double multipliersWQ21[]={1, 3, 20, 50,
                        200, 1, 1, 1,
                        1, 1, 1, 1};
  double multipliersxy1[]={1, 5, 1, 1,
                        1, 2, 1, 1,
                        1, 1, 1, 1};

  //vector<double> multipliers = GetScales(histAndOpts,false,5,0.75,false);
  GridCanvas* gc=plotXAxis1D(histAndOpts, "Muon Energy (GeV)", "Ehad", doMultipliers ? multipliers1 : NULL);
  //GridCanvas* gc=plotXAxis1D(histAndOpts, "Muon Energy (GeV)", "Ehad", doMultipliers ? &multipliers[0] : NULL);
  //GridCanvas* gc=plotXAxis1D(histAndOpts, "x", "Q2", doMultipliers ? multipliersxQ21 : NULL);
  //GridCanvas* gc=plotXAxis1D(histAndOpts, "W", "Q2", doMultipliers ? multipliersWQ21 : NULL);
  //GridCanvas* gc=plotXAxis1D(histAndOpts, "Bjorken Variable x", "y", doMultipliers ? multipliersxy1 : NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  //gc->Remax();
  gc->SetYLimits(0, 2.0);
  //gc->SetYLimits(0, 90);
  //gc->SetYLimits(0, 15);
  //gc->SetYLimits(0, 200);
  gc->SetYTitle("Event Rate(x10^{-3}) per GeV^{2}");
  //gc->SetYTitle("Event Rate");
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
  leg->AddEntry(TinTrans, "Trans in Trans", "l");
  leg->AddEntry(CinTrans, "Contin in Trans", "l");
  leg->AddEntry(SinTrans, "SinTrans", "l");
  leg->AddEntry(TransUS, "US in Trans", "l");
  leg->AddEntry(TransDS, "DS in Trans", "l");
  leg->AddEntry(" ","Iron of target 2");
  //leg->AddEntry(" ","Lead of Target 3");
  //leg->AddEntry(" ","Lead of Target 3");
  leg->Draw("SAME");
    gc->Print(doMultipliers ? "nu-2d-evtrate-model-Ehad-multiplier.png" : "nu-2d-evtrate-model-Ehad.png");
    //gc->Print(doMultipliers ? "nu-2d-evtrate-model-Q2-multiplier.png" : "nu-2d-evtrate-model-x.png");
    //gc->Print(doMultipliers ? "nu-2d-evtrate-model-Q2-multiplier.png" : "nu-2d-evtrate-model-W.png");
    //gc->Print(doMultipliers ? "nu-2d-evtrate-model-y-multiplier.png" : "nu-2d-evtrate-model-y.png");
  


  // ------------------------------------------------------------------------------
  //
  // Now make Emu in bins of Ehad. It's all the same

  // Values to multiply each bin by to get them on a similar range
  double  multipliers4[]={4, 1, 1, 1, 1.5, 3, 5.5, 8, 12, 20,
			  40, 100, 1};
  double  multipliersxQ22[]={3, 1, 1, 2, 10,
			  100, 1, 1, 1, 1,
			  1, 1, 1};
  double  multipliersWQ22[]={1, 2, 7, 50, 500,
			  1, 1, 1, 1, 1,
			  1, 1, 1};

  double  multipliersxy2[]={4, 1.5, 1, 1, 3,
			  20, 1, 1, 1, 1,
			  1, 1, 1};
  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  //vector<double> multipliers2 = GetScales(histAndOpts,true,5,0.75,false);
  GridCanvas* gc2=plotYAxis1D(histAndOpts, "Hadronic Energy (GeV)", "Emu", doMultipliers ? multipliers4 : NULL);
  //GridCanvas* gc2=plotYAxis1D(histAndOpts, "Hadronic Energy (GeV)", "Emu", doMultipliers ? &multipliers2[0] : NULL);
  //GridCanvas* gc2=plotYAxis1D(histAndOpts, "Q2 (GeV^{2})", "x", doMultipliers ? multipliersxQ22 : NULL);
  //GridCanvas* gc2=plotYAxis1D(histAndOpts, "Q2 (GeV^{2})", "W", doMultipliers ? multipliersWQ22 : NULL);
  //GridCanvas* gc2=plotYAxis1D(histAndOpts, "Inelasticity y", "x", doMultipliers ? multipliersxy2 : NULL);
  gc2->SetYLimits(0, 0.6);
  //gc2->SetYLimits(0, 90);
  //gc2->SetYLimits(0, 15);
  //gc2->SetYLimits(0, 200);
  //gc2->Remax();
  gc2->SetYTitle("Event Rate(x10^{-3}) per GeV^{2}");
  //gc2->SetYTitle("Event Rate");
  gc2->Modified();
    

  gc2->Print(doMultipliers ? "nu-2d-evtrate-model-Emu-multiplier.png" : "nu-2d-evtrate-model-Emu.png");
  //gc2->Print(doMultipliers ? "nu-2d-evtrate-model-x-multiplier.png" : "nu-2d-evtrate-model-Q2.png");
  //gc2->Print(doMultipliers ? "nu-2d-evtrate-model-W-multiplier.png" : "nu-2d-evtrate-model-Q2.png");
  //gc2->Print(doMultipliers ? "nu-2d-evtrate-model-x-multiplier.png" : "nu-2d-evtrate-model-x.png");

}

int main(int argc, char* argv[])
{

  string location = argv[1];
  //multipliers
    makePlots(true);
  //makePlots(true,true,true,location);
  //makePlots(true,false,false,location);
  //standard
/*  makePlots(false,true,false,location);
  makePlots(false,true,true,location);
  makePlots(false,false,false,location);
*/
  return 0;
}
