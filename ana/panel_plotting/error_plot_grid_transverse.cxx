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

#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

#include <iostream>

using namespace std;
using namespace PlotUtils;

// =====================================================================
// axis==1: x axis. axis==2: y axis

// Take a vector of 1D histograms which all have the same binning and
// stack them together into a 2D histogram. The axis argument says
// which axis to stack them on
TH2* concatenateHists(vector<TH1*>& hists1D, int axis, TH2D* temp)
{
  assert(hists1D.size());


  TH2* ret=0;
  if(axis==1){
    ret=new TH2D(uniq(), TString::Format(";%s", hists1D[0]->GetXaxis()->GetTitle()),
                 hists1D[0]->GetXaxis()->GetNbins(), hists1D[0]->GetXaxis()->GetXbins()->GetArray(),
                 temp->GetYaxis()->GetNbins(), temp->GetYaxis()->GetXbins()->GetArray());
  }
  else{
    ret=new TH2D(uniq(), TString::Format(";;%s", hists1D[0]->GetXaxis()->GetTitle()),
                 temp->GetXaxis()->GetNbins(), temp->GetXaxis()->GetXbins()->GetArray(),
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

  return ret;
}

// =====================================================================
vector<std::pair<TH2*, const char*> > getSystHistsAndOpts(MnvH2D* data, bool pt, TLegend *&leg)
{
  MnvPlotter plotter;
  plotter.ApplyStyle(kCCQENuStyle);
  
  // For each bin in the other variable, make a vector of the
  // systematic histograms
  const int nBins=pt ? data->GetXaxis()->GetNbins() : data->GetYaxis()->GetNbins();
  vector<vector<TH1*> > histsPT;
  histsPT.resize(nBins);

  // Get MnvPlotter to plot all the histograms, and slurp them into histsPT
  for(int i=0; i<nBins; ++i){
    // First plot the histograms in the dummy canvas...
    TCanvas c;
    MnvH1D* proj=pt ? data->ProjectionY(uniq(), i+1, i+1) : data->ProjectionX(uniq(), i+1, i+1);
    plotter.DrawErrorSummary(proj, "TR", true, true, -1, false, "", true);
    leg = new TLegend(*getPadLegend(&c));
    leg->SetName("MyLegend");
    //    c.Print("can.png");
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
    for(unsigned int iBin=0; iBin<histsPT.size(); ++iBin){  
      histsPT_transpose[iSyst].push_back(histsPT[iBin][iSyst]);
    }
  }
  cout << "Loop over syst done" << endl;
  vector<std::pair<TH2*, const char*> > histsPT2D;
  // TODO: Figure out why the last systematic is crashing
  for(int iSyst=0; iSyst<histsPT_transpose.size(); ++iSyst){
    TH2* h2d=concatenateHists(histsPT_transpose[iSyst], pt ? 2 : 1,data);
    // We want to draw all of these histograms as graphs, and exclude
    // the zero bins, to get rid of ROOT artifacts. The "graph0" draw
    // option does that (and I made it safe to pass all graphs)
    histsPT2D.push_back(std::make_pair(h2d, "graph0 l"));
  }
  cout << "Made the concat is now done\t" << histsPT2D.size() << endl;
  return histsPT2D;
}

// =====================================================================
void makePlots(bool pt, string location, string variablex, string variabley="", int projx_column=-1, int projx_row=-1, int projx_pixelx=-1, int projx_pixely=-1, int projy_column=-1, int projy_row=-1, int projy_pixelx=-1, int projy_pixely=-1, string sample="")
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
  cout << location << "\t" << variablex << "\t" << variabley << "\t" << endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  string variable = "";
  string hvariable = "";
  string fname = "";
  if(variabley!=""){
    variable = variablex+"_"+variabley;
    hvariable = variablex+"_"+variabley;
    fname =Form("%s_CV/CrossSection_%s_CombinedPlaylists.root",location.c_str(),variable.c_str()) ;
  }
  else{
    variable = variablex;
    hvariable = variablex+"_ptmu";
    variabley="ptmu";
    fname = Form("%s",location.c_str());
  }
  string histname = Form("h_%s_data_nobck_unfold_effcor_cross_section",hvariable.c_str());
  if(variabley!="ptmu"){
    if(sample=="CrossSection")histname = Form("h_%s_data_nobck_unfold_effcor_cross_section",hvariable.c_str());
    else if(sample=="Efficiency") histname = Form("h_%s_data_nobck_unfold_effcor",hvariable.c_str());
    else if(sample=="Unfolding") histname = Form("h_%s_data_nobck_unfold",hvariable.c_str());
    else if(sample=="BkgSub") histname = Form("h_%s_data_nobck",hvariable.c_str());
    else exit(1);
  }

  TFile f(fname.c_str());
  MnvH2D* data=(MnvH2D*)f.Get(histname.c_str());
  cout << "\t" << data << endl;
  TLegend *leg = NULL;
  vector<std::pair<TH2*, const char*> > histsPT2D=getSystHistsAndOpts(data,  pt, *& leg);
  leg->SetLineColor(kWhite);
  leg->SetFillColor(kWhite);
  cout << "PLOT PLOT PLOT\t"  << data<< endl;
  if(pt){
    cout << "Draw PT" << endl;
    GridCanvas* gcPT;
    if(projx_column!=-1) gcPT=plotYAxis1D(histsPT2D,variabley.c_str(),variablex.c_str(),projy_column,projy_row,projy_pixelx,projy_pixely);
    else gcPT=plotYAxis1D(histsPT2D,variabley.c_str(),variablex.c_str());
    cout << "HERE" << endl;
    gcPT->SetYTitle("Fractional uncertainty");
    gcPT->SetYLimits(0, 0.29);
    if(projy_column!=-1){
      leg->SetX1(0.02);
      leg->SetY1(0.92);
      leg->SetX2(0.98);
      leg->SetY2(1.0);
      leg->SetNColumns(3);
    }
    else{
      leg->SetX1(0.58);
      leg->SetY1(0.05);
      leg->SetX2(0.95);
      leg->SetY2(0.32);
    }      
    leg->Draw("SAME");
    gcPT->Modified();
    gcPT->Print("errors-pt.eps");
    gcPT->Print("errors-pt.png");
    gcPT->Print("errors-pt.C");
  }
  else{
    GridCanvas* gcPT;
    if(projy_column!=-1) gcPT=plotXAxis1D(histsPT2D,variablex.c_str(),variabley.c_str(),projx_column,projx_row,projx_pixelx,projx_pixely);
    else gcPT=plotXAxis1D(histsPT2D,variablex.c_str(),variabley.c_str());
    gcPT->SetYTitle("Fractional uncertainty");
    gcPT->SetYLimits(0,0.29);
    if(projy_column!=-1){
      leg->SetX1(0.02);
      leg->SetY1(0.92);
      leg->SetX2(0.98);
      leg->SetY2(1.0);
      leg->SetNColumns(3);
    }
    else{
      leg->SetX1(0.58);
      leg->SetY1(0.05);
      leg->SetX2(0.95);
      leg->SetY2(0.32);
    }      
    leg->Draw("SAME");
    gcPT->Modified();
    gcPT->Print("errors-pz.eps");
    gcPT->Print("errors-pz.png");
    gcPT->Print("errors-pz.C");
  }
}

int main(int argc, char* argv[])
{
  if(argc==3){
    cout << argv[1] << "\t" << argv[2] << endl;
    makePlots(true,argv[1],argv[2]);
    makePlots(false,argv[1],argv[2]);
  }
  else if(argc>3){
    makePlots(false,argv[1],argv[2],argv[3],atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9]),atoi(argv[10]),atoi(argv[11]),argv[12]);
    makePlots(true,argv[1],argv[2],argv[3],atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9]),atoi(argv[10]),atoi(argv[11]),argv[12]);
  }
  else return 1;
  return 0;
}
