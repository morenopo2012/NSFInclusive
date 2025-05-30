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
    double scale = 3.5/maxval;
    if(scale>1){
      int tmpscale = floor(scale*10);
      scale = tmpscale/10.0;
    }
    else{
      int tmpscale = ceil(scale*10);
      scale = tmpscale/10.0;
    }
    cout << scale << endl;
    tmpvect.push_back(scale);
  }
    return tmpvect;
}


void makePlots(bool doMultipliers,bool doGenies,string location,string varx, string vary)
{
//  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  string var = varx+"_"+vary;
  string varspecfilename = var;
  string varspec = var;
  if(vary=="ptmu") varspecfilename=varx;
  cout << varspecfilename << endl;
  string filename = Form("%s_CV_Elastic_FSI//selection_Signal/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_CombinedPlaylists.root",location.c_str());//Final result
  //  if(vary!="ptmu")filename = Form("%s_CV/CrossSection_%s_CombinedPlaylists.root",location.c_str(),varspecfilename.c_str());//Final result
  TFile f1(filename.c_str());
  //  TFile f2(Form("%s_default/CrossSection_CombinedPlaylists.root_%s.root",location.c_str(),varspecfilename.c_str()));//Default GENIE
  //  TFile f3(Form("%s_pion_rpa/CrossSection_CombinedPlaylists.root_%s.root",location.c_str(),varspecfilename.c_str()));//GENIE+2p2h+RPA (aka no tune)
  MnvH2D* dataMnv=(MnvH2D*)f1.Get(Form("h_%s_data",varspec.c_str()));
  MnvH2D* mcMnv=(MnvH2D*)f1.Get(Form("h_%s_mc",varspec.c_str()));

  MnvH2D* mcMnv_qelike_qe = (MnvH2D*)f1.Get(Form("h_%s_qelike_qe",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_res = (MnvH2D*)f1.Get(Form("h_%s_qelike_res",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_dis = (MnvH2D*)f1.Get(Form("h_%s_qelike_dis",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_2p2h = (MnvH2D*)f1.Get(Form("h_%s_qelike_2p2h",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelikenot   = (MnvH2D*)f1.Get(Form("h_%s_qelikenot",varspec.c_str()));//Get from N track



  dataMnv->Scale(1e-5, "width");
  mcMnv->Scale(1e-5, "width");
  mcMnv_qelike_qe->Scale(1e-5,"width");
  mcMnv_qelike_res->Scale(1e-5,"width");
  mcMnv_qelike_dis->Scale(1e-5,"width");
  mcMnv_qelike_2p2h->Scale(1e-5,"width");
  mcMnv_qelikenot->Scale(1e-5,"width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());

  TH2* mc_qelike_qe = new TH2D(mcMnv_qelike_qe->GetCVHistoWithStatError());
  TH2* mc_qelike_res = new TH2D(mcMnv_qelike_res->GetCVHistoWithStatError());
  TH2* mc_qelike_dis = new TH2D(mcMnv_qelike_dis->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h = new TH2D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());
  TH2* mc_qelikenot  = new TH2D(mcMnv_qelikenot->GetCVHistoWithStatError());

  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);
  mc_qelike_qe->SetLineColor(mycolors[3]);
  mc_qelike_res->SetLineColor(mycolors[4]);
  mc_qelike_dis->SetLineColor(mycolors[5]);
  mc_qelike_2p2h->SetLineColor(mycolors[16]);
  mc_qelikenot->SetLineColor(mycolors[8]);

  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.5);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);

  dataStat->SetLineColor(kBlack);


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
  //  histAndOpts.push_back(std::make_pair(dataStat, "histpe1"));
  histAndOpts.push_back(std::make_pair(mc,       "hist l"));
  histAndOpts.push_back(std::make_pair(mc_qelike_qe,       "hist l"));
  histAndOpts.push_back(std::make_pair(mc_qelike_res,       "hist l"));
  histAndOpts.push_back(std::make_pair(mc_qelike_dis,       "hist l"));
  histAndOpts.push_back(std::make_pair(mc_qelike_2p2h,       "hist l"));
  histAndOpts.push_back(std::make_pair(mc_qelikenot,       "hist l"));
  //  histAndOpts.push_back(std::make_pair(data,     "histpe1"));



  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  vector<double> multi_1 = GetScales(histAndOpts, false);

  GridCanvas* gc=NULL;
  gc= plotYAxis1D(histAndOpts, vary.c_str() ,varx.c_str() ,doMultipliers ? &multi_1[0] : NULL);
  

  // Set the y range manually. Can also use gc->Remax() to guess automatically
  gc->SetYLimits(0, 4.99);
  
  //Label thy axis
  gc->SetYTitle("Event Rate(*1e-5) per GeV^{2}");
  gc->SetXTitle("Reconstructed Struck Nucleon Momentum (GeV)");
    
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  //  TLegend* leg=new TLegend(0.17, 0.7, 0.31, 0.9);
  TLegend* leg=new TLegend(0.4, 0.1, 0.9, 0.3);
  leg->SetNColumns(2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MnvGENIEv1.0.1", "l");
  leg->AddEntry(mc_qelike_qe,"QE","l");
  leg->AddEntry(mc_qelike_res,"Resonant","l");
  leg->AddEntry(mc_qelike_dis,"DIS","l");
  leg->AddEntry(mc_qelike_2p2h,"2p2h","l");
  leg->AddEntry(mc_qelikenot,"Background","l");
  leg->Draw("SAME");

  gc->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projy-%s-multiplier.eps",var.c_str()) : Form("nu-2d-evtrate-comps-projy-%s.eps",var.c_str()));
  gc->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projy-%s-multiplier.png",var.c_str()) : Form("nu-2d-evtrate-comps-projy-%s.png",var.c_str()));
  gc->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projy-%s-multiplier.C",var.c_str()) : Form("nu-2d-evtrate-comps-projy-%s.C",var.c_str()));


  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  vector<double> multi_2 = GetScales(histAndOpts,true);
  


  // Values to multiply each bin by to get them on a similar range


  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=NULL;
  gc2 =plotXAxis1D(histAndOpts, varx.c_str() , vary.c_str() ,doMultipliers ? &multi_2[0] : NULL);
  gc2->SetYLimits(0, 4.99);
  gc2->SetYTitle("Event rate(x1e-5) per GeV^{2}");
  gc2->SetXTitle("Reconstructed Struck Nucleon Momentum (GeV)");
  leg->Draw("SAME");
  gc2->Modified();
  
  gc2->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projx-%s-multiplier.eps",var.c_str()) : Form("nu-2d-evtrate-comps-projx-%s.eps",var.c_str()));
  gc2->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projx-%s-multiplier.png",var.c_str()) : Form("nu-2d-evtrate-comps-projx-%s.png",var.c_str()));
  gc2->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projx-%s-multiplier.C",var.c_str()) : Form("nu-2d-evtrate-comps-projx-%s.C",var.c_str()));


}

int main(int argc, char* argv[])
{
  // makePlots(true,true,argv[1],argv[2]);
  // makePlots(true,false,argv[1],argv[2]);
  // makePlots(false,true,argv[1],argv[2]);
  //  makePlots(false,false,argv[1],argv[2],argv[3]);
  makePlots(true,false,argv[1],argv[2],argv[3]);

  return 0;
}
