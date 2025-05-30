//  #include "../util/plot/myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvPlotter.h"

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

//#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

#include <iostream>
#include <algorithm>
using namespace std;
using namespace PlotUtils;

// =====================================================================
// axis==1: x axis. axis==2: y axis

// Take a vector of 1D histograms which all have the same binning and
// stack them together into a 2D histogram. The axis argument says
// which axis to stack them on
TH2* concatenateHists(vector<TH1*>& hists1D, int axis)
{
  assert(hists1D.size());

 cout << "concatenateHists with " << hists1D.size() << " hists, axis=" << axis << endl;
  const int nyBins=7;
  const double yBins[nyBins+1]={1.5, 3., 5., 8., 11., 15., 20., 25. };

  const int nxBins=11;
  const double xBins[nxBins+1]={2.0, 3.5, 4.75, 5.5, 7.5, 10., 13., 16., 20., 25., 35., 50 };

  TH2* ret=0;
  if(axis==1){
    ret=new TH2D(uniq(), TString::Format(";%s", hists1D[0]->GetXaxis()->GetTitle()),
                 hists1D[0]->GetXaxis()->GetNbins(), hists1D[0]->GetXaxis()->GetXbins()->GetArray(),
                 nyBins, yBins);
  }
  else{
    ret=new TH2D(uniq(), TString::Format(";;%s", hists1D[0]->GetXaxis()->GetTitle()),
                 nxBins, xBins,
                 hists1D[0]->GetXaxis()->GetNbins(), hists1D[0]->GetXaxis()->GetXbins()->GetArray());
  }

  ret->SetLineColor(hists1D[0]->GetLineColor());
  ret->SetLineStyle(hists1D[0]->GetLineStyle());

  for(unsigned int iHist=0; iHist<hists1D.size(); ++iHist){
    for(int j=0; j<hists1D[0]->GetXaxis()->GetNbins()+1; ++j){
      int ixBin=axis==1 ? j       : iHist+1;
      int iyBin=axis==1 ? iHist+1 : j;
      double content=hists1D[iHist]->GetBinContent(j);
      ret->SetBinContent(ixBin, iyBin, content);
    }
  }
  cout << "This is my concat hist address " << ret << endl;
  return ret;
}

// =====================================================================
vector<std::pair<TH2*, const char*> > getSystHistsAndOpts(MnvH2D* data, bool pt, string group = "")
{
  MnvPlotter plotter;
  plotter.ApplyStyle(kCCQENuInclusiveStyle);
  vector<string> vertnames = data->GetVertErrorBandNames();
  vector<string> latnames = data->GetLatErrorBandNames();
  // For each bin in the other variable, make a vector of the
  // systematic histograms
  const int nBins=pt ? 16 : 15;
  vector<vector<TH1*> > histsPT;
  histsPT.resize(nBins);

  // Get MnvPlotter to plot all the histograms, and slurp them into histsPT
  for(int i=0; i<nBins; ++i){
    // First plot the histograms in the dummy canvas...
    TCanvas c;
    MnvH1D* proj=pt ? data->ProjectionY(uniq(), i+1, i+1) : data->ProjectionX(uniq(), i+1, i+1);
    if(group=="") plotter.DrawErrorSummary(proj, "N", true, true, -1, false,"", true);
    else{
      if(std::count(vertnames.begin(),vertnames.end(),group)){
	TH1D err = proj->GetVertErrorBand(group)->GetErrorBand(true);
	err.DrawClone();
      }
      else{
	TH1D err = proj->GetLatErrorBand(group)->GetErrorBand(true);
	err.DrawClone();
      }
    }
    std::vector<TH1*> padHists=getPadHists(&c);
    histsPT[i]=padHists;
  }

  // concatenateHists wants a vector of hists for each of the bins of
  // a given systematic. But histsPT is the other way round (the inner
  // vector loops over systematics).  So we have this fiddly loop to
  // make a transposed version of the vector-of-vector

  // It would have been easier to just pass the original
  // vector-of-vector into concatenateHists, and tell it which
  // systematic we wanted, but I've written and debugged this now, so
  // not changing it

  //  First index is systematic, second
  // index is bin
  vector<vector<TH1*> > histsPT_transpose;
  int nSyst=histsPT[0].size();
  cout << "There are " << nSyst << " systematics" << endl;
  histsPT_transpose.resize(nSyst);

  for(int iSyst=0; iSyst<nSyst; ++iSyst){
    cout << histsPT[0][iSyst]->GetName() << endl;
    for(unsigned int iBin=0; iBin<histsPT.size(); ++iBin){  
      histsPT_transpose[iSyst].push_back(histsPT[iBin][iSyst]);
    }
  }

  vector<std::pair<TH2*, const char*> > histsPT2D;
  // TODO: Figure out why the last systematic is crashing
  for(int iSyst=0; iSyst<histsPT_transpose.size(); ++iSyst){
    cout << iSyst << endl;
    TH2* h2d=concatenateHists(histsPT_transpose[iSyst], pt ? 2 : 1);
    // We want to draw all of these histograms as graphs, and exclude
    // the zero bins, to get rid of ROOT artifacts. The "graph0" draw
    // option does that (and I made it safe to pass all graphs)
    histsPT2D.push_back(std::make_pair(h2d, "graph0 l"));
  }
  cout << "Done getSystHistsAndOpts " << endl;
  return histsPT2D;
}

// =====================================================================
void makePlots(bool pt, string location, bool drawGroups, int areanorm)
{
  // This turns out to be complicated (well, I could have made it less
  // bad, but this way is complicated and generalizable):
  //
  // MnvPlotter knows how to make the histograms we need, including
  // fancy stuff like grouping systematics together and sticking to
  // colour schemes. But it doesn't know about GridCanvas, and it
  // goes off and draws the histograms in its own way
  //
  // So we're going to take projections, and ask MnvPlotter to plot
  // the projection in a dummy canvas. Then we'll grab all the
  // histograms from that canvas and hold onto them.
  //
  // We could just grab all those 1D histograms and plot them straight
  // into the appropriate panel of our GridCanvas, but I want to reuse
  // the plotpz1D and plotpT1D functions in plot.h, because they know
  // how to do fanciness like squashing the tail in pz. Those
  // functions take a vector of 2D histograms, so we have to take our
  // 1D histograms from MnvPlotter, and stack them back together into
  // 2D histograms. That's what the concatenateHists function does

  //ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  //gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  //TFile f(Form("%s_CV/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));
  //TFile f1("histsFullStatDIS/Hists_EventSelection_t1_z26_Nu_v1_.root");//1track
  TFile f1("/minerva/data/users/zdar/MADHists/histsFullStatDIS/histsFullStatDIS/Hists_EventSelection_t2_z26_Nu_v1_.root");//1track
  
  MnvH2D* data=(MnvH2D*)f1.Get("selected_mc_reco2d_Emu_Ehad");

  MnvH2D* zerobin = (MnvH2D*)data->Clone("zero");
  zerobin->Reset();
  zerobin->ClearAllErrorBands();
  for(int i=0;i<zerobin->GetNbinsX();i++){
    for(int j=0;j<zerobin->GetNbinsY();j++){
      if(i==5 && j==12){
	zerobin->SetBinContent(i+1,j+1,0.0);
	zerobin->SetBinError(i+1,j+1,0.0);
      }
      else{
	zerobin->SetBinContent(i+1,j+1,1.0);
	zerobin->SetBinError(i+1,j+1,0.0);
      }
      
    }
  }
  zerobin->AddMissingErrorBandsAndFillWithCV(*data);
  zerobin->SaveAs("test.root");
  //  data->Multiply(data,zerobin);
  if(areanorm==1) data = new MnvH2D(data->GetAreaNormalizedCopy());
  vector<std::pair<TH2*, const char*> > histsPT2D=getSystHistsAndOpts(data,  pt);
  vector<string> names = data->GetErrorBandNames();
  cout << "Start Sysnames " << endl;
  if(!drawGroups){
    vector<string> Sysnames;
    Sysnames.push_back("Total Uncertainty");
    Sysnames.push_back("Statisical");
    Sysnames.push_back("Flux");
    //Sysnames.push_back("Hadrons");
    Sysnames.push_back("Models");
    Sysnames.push_back("Muon Reconstruction");
    Sysnames.push_back("MINOS Efficiency");
    //Sysnames.push_back("Normalization");


    
    
    
    TLegend *leg = new TLegend(0.6, 0.05, 0.95, 0.32);
    leg->SetNColumns(2);
    leg->SetLineColor(kWhite);
    leg->SetFillColor(kWhite);
    for(int i=0;i<histsPT2D.size();i++){
      cout << histsPT2D[i].first->GetName() << endl;
      leg->AddEntry(histsPT2D[i].first,Sysnames[i].c_str(),"l");
    }
    cout << "I'm here" << endl;
    if(pt){
      //GridCanvas* gcPT=plotXAxis1DRebinPz(histsPT2D, "P_{||} (GeV)", "P_{t}", 4,4,800,500,NULL);
      GridCanvas* gcPT=plotXAxis1DRebinPz(histsPT2D, "Muon Energy (GeV)", "Ehad", 4,4,800,500,NULL);
      gcPT->SetYTitle("Fractional uncertainty");
      gcPT->Remax();

      leg->Draw("SAME");
      gcPT->Modified();
      //gcPT->Print(Form("errors-pt_areanorm_%d.eps",areanorm));
      gcPT->Print(Form("errors-Emu_areanorm_%d.png",areanorm));
      //gcPT->Print(Form("errors-pt_areanorm_%d.C",areanorm));
    }
    else{
      //GridCanvas* gcPT=plotYAxis1D(histsPT2D, "P_{t} (GeV)", "P_{||}",4,4,800,500,NULL);
      GridCanvas* gcPT=plotYAxis1D(histsPT2D, "Hadronic Energy (GeV)", "Emu",4,4,800,500,NULL);
      gcPT->SetYTitle("Fractional uncertainty");
      gcPT->Remax();
      gcPT->Modified();
      //gcPT->Print(Form("errors-pz_areanorm_%d.eps",areanorm));
      gcPT->Print(Form("errors-Ehad_areanorm_%d.png",areanorm));
      //gcPT->Print(Form("errors-pz_areanorm_%d.C",areanorm));
    }
  }
  else{
    for(int n=0;n<names.size();n++){
      string group = names[n];
      vector<std::pair<TH2*, const char*> > histsPT2D=getSystHistsAndOpts(data,  pt, group);   
      TLegend *leg = new TLegend(0.78, 0.05, 0.95, 0.32);
      leg->SetLineColor(kWhite);
      leg->SetFillColor(kWhite);
      for(int j=0;j<histsPT2D.size();j++)leg->AddEntry(histsPT2D[j].first,group.c_str(),"l");
      if(pt){
	//GridCanvas* gcPT=plotXAxis1DRebinPz(histsPT2D, "P_{||} (GeV)", "P_{t}", 4,4,800,500,NULL);
	GridCanvas* gcPT=plotXAxis1DRebinPz(histsPT2D, "Muon Energy (GeV)", "Ehad", 4,4,800,500,NULL);
	gcPT->SetYTitle("Fractional uncertainty");
        gcPT->Remax();
 	leg->Draw("SAME");
	gcPT->Modified();
	//gcPT->Print(Form("errors-%s-pt_areanorm_%d.eps",group.c_str(),areanorm));
	gcPT->Print(Form("errors-%s-Emu_areanorm_%d.png",group.c_str(),areanorm));
	//gcPT->Print(Form("errors-%s-pt_areanorm_%d.C",group.c_str(),areanorm));
      }
      else{
	//GridCanvas* gcPT=plotYAxis1D(histsPT2D, "P_{t} (GeV)", "P_{||}",4,4,800,500,NULL);
	GridCanvas* gcPT=plotYAxis1D(histsPT2D, "Hadronic Energy (GeV)", "Emu",4,4,800,500,NULL);
	gcPT->SetYTitle("Fractional uncertainty");
	//gcPT->SetYLimits(0, 0.25);
        gcPT->Remax();
	//gcPT->SetYLimits(0, 0.29);
	gcPT->Modified();
	//gcPT->Print(Form("errors-%s-pz_areanorm_%d.eps",group.c_str(),areanorm));
	gcPT->Print(Form("errors-%s-Ehad_areanorm_%d.png",group.c_str(),areanorm));
	//gcPT->Print(Form("errors-%s-pz_areanorm_%d.C",group.c_str(),areanorm));
      }
    }
  }
}

int main(int argc, char* argv[])
{
  for(int i=0;i<2;i++){
    makePlots(true,argv[1],false,i);
    makePlots(false,argv[1],false,i);
    
    makePlots(true,argv[1],true,i);
    makePlots(false,argv[1],true,i);
  }
  return 0;
}
