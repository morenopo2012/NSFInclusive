//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/HyperDimLinearizer.h"//THIS HAS TO CHANGE TO BE INCLUDED IN THE MAKE FILE EVENTUALLY.
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
#include "localColor.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"
#include <iostream>
#include <fstream>
#include <string>

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
    double scale = 75000/maxval;
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



void makePlots(bool doMultipliers, string location,string sample, bool doRatio, string varx, string vary)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  MnvPlotter *plotter = new MnvPlotter();

  string scaleform = "Signal";
  //Old
  /*
  if(sample=="MichelSideBand") scaleform = "_michel";
  if(sample=="BlobSideBand") scaleform = "_blobs";
  if(sample=="MicBlobSideBand") scaleform = "_micblob";
  */
  if(sample=="MichelSideBand") scaleform = "Michel";
  if(sample=="BlobSideBand") scaleform = "Blob";
  if(sample=="MicBlobSideBand") scaleform = "Micblob";

  string sbstring =   varx+"_"+vary;

  cout << scaleform << "\t" <<Form("%s/selection_%s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-%s_CombinedPlaylists.root",location.c_str(),sample.c_str(),sample.c_str()) << endl;
  TFile f1(Form("%s/selection_%s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-%s_CombinedPlaylists.root",location.c_str(),sample.c_str(),sample.c_str()));//sample 1 trk
  TFile f2(Form("%s/SideBandFit_%s_CombinedPlaylists.root",location.c_str(),sbstring.c_str()));//sample 1 trk
  string histstring = varx+"_"+vary;
  
  cout<<histstring << endl;
  MnvH2D* dataMnv=(MnvH2D*)f1.Get(Form("h_%s_data",histstring.c_str()));
  MnvH2D* mcqenotMnv=(MnvH2D*)f1.Get(Form("h_%s_qelikenot",histstring.c_str()));
  MnvH2D* mcqenotscaledMnv=mcqenotMnv->Clone("Constrainedqelikenot");

  MnvH2D* mcqeMnv=(MnvH2D*)f1.Get(Form("h_%s_qelike",histstring.c_str()));
  MnvH2D* mcqescaledMnv=mcqeMnv->Clone("Constrainedqelike");

  string vary_mod;
  if(vary=="ptmu") vary_mod="ptmu";
  else vary_mod=vary;
  string bkgstring=varx+vary_mod;
  //This is the old ptmu constraint method
  //  MnvH2D* scalefactor = (MnvH2D*)f2.Get(Form("h_weights_2track_%sbins_qelikenot%s",bkgstring.c_str(),scaleform.c_str()));
  MnvH2D* scalefactor = (MnvH2D*)f2.Get(Form("%s_qelikenot_scaled",scaleform.c_str()));
  mcqenotscaledMnv->Multiply(mcqenotscaledMnv,scalefactor);
  //UnScaled QELike
  MnvH2D* mcMnv = mcqenotMnv->Clone("mcMnv_unconstraint");
  MnvH2D* mcscaledMnv =mcqenotscaledMnv->Clone("mcMnv_constraint");
  mcMnv->Add(mcqeMnv);
  mcscaledMnv->Add(mcqeMnv);




  int ndf = 0;
  double chi2_unconstrained = plotter->Chi2DataMC(dataMnv,mcMnv,ndf,1,false,false,true,NULL);
  double chi2_constrained = plotter->Chi2DataMC(dataMnv,mcscaledMnv,ndf,1,false,false,true,NULL);

  TH2* dataresults = new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mcresults = new TH2D(mcMnv->GetCVHistoWithError());
  TH2* mcresultsscaled = new TH2D(mcscaledMnv->GetCVHistoWithError());
  
  
  vector<int> mycolors = getColors(2);
  mcresults->SetLineColor(kRed);
  mcresults->SetLineWidth(2);
  mcresults->SetFillColor(kRed);
  mcresults->SetFillStyle(3005);

  mcresultsscaled->SetLineColor(kBlue);
  mcresultsscaled->SetLineWidth(2);
  mcresultsscaled->SetLineStyle(2);
  mcresultsscaled->SetFillColor(kBlue);
  mcresultsscaled->SetFillStyle(3004);

  
  dataresults->SetMarkerStyle(kFullCircle);
  dataresults->SetMarkerSize(0.7);
  dataresults->SetLineColor(kBlack);
  dataresults->SetLineWidth(2);


  dataresults->Scale(1,"width");
  mcresults->Scale(1,"width");
  mcresultsscaled->Scale(1,"width");
  


  if(doRatio){
    dataresults->Divide(mcresults);
    mcresultsscaled->Divide(mcresults);
    mcresults->Divide(mcresults);
  }
  //  dataresults->GetXaxis()->SetRangeUser(0,1000);
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
  histAndOpts.push_back(std::make_pair(dataresults, "histpe"));
  histAndOpts.push_back(std::make_pair(mcresults,       "histle3"));
  histAndOpts.push_back(std::make_pair(mcresultsscaled,       "histle3"));
  histAndOpts.push_back(std::make_pair(dataresults, "histpe"));
  cout << "Mults" << endl;
  

  vector<double> multi_1 = GetScales(histAndOpts, true);
  
    
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0, 0, 1, 1);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(dataresults, "MINERvA data", "lpe");
  leg->AddEntry(mcresults, "Minerva Tune v1.0.1", "l");
  leg->AddEntry(mcresultsscaled, "Minerva Tune v1.0.1 constrained", "l");
  
  
  TLatex mytex;
  mytex.SetTextSize(0.05);
  string mystring = Form("Chi2 before: %0.2f ; Chi2 after: %0.2f ; ndf = %d",chi2_unconstrained,chi2_constrained,ndf);
  
  
  GridCanvas* gc1=NULL;
  gc1=plotXAxis1D(histAndOpts, varx,vary,doMultipliers ? &multi_1[0] : NULL);
  
  // Set the y range manually. Can also use gc1->Remax() to guess automatically
  if(doRatio)gc1->SetYLimits(0,1.99);
  else gc1->SetYLimits(0, 100000);
  //gc1->SetYLimits(0, 1e-39);
  if(doRatio)gc1->SetYTitle("Ratio to MnvGENIEv1");
  else gc1->SetYTitle(Form("d^{2}N/d%sd%s",varx.c_str(),vary.c_str()));
  mytex.DrawLatex(0.04,0.96,mystring.c_str());
  gc1->Modified();
  
  if(doRatio){
    gc1->Print(doMultipliers ? Form("nu-tki-ptmuconstraint_%s-comps-pt-multiplier_ratio-2track_%s_%s.eps",sample.c_str(),varx.c_str(),vary.c_str()) : Form("nu-tki-ptmuconstraint-%s-comps-pt_ratio-2track_%s_%s.eps",sample.c_str(),varx.c_str(),vary.c_str()));
    gc1->Print(doMultipliers ? Form("nu-tki-ptmuconstraint-%s-comps-pt-multiplier_ratio-2track_%s_%s.png",sample.c_str(),varx.c_str(),vary.c_str()) : Form("nu-tki-ptmuconstraint-%s-comps-pt_ratio-2track_%s_%s.png",sample.c_str(),varx.c_str(),vary.c_str()));
    gc1->Print(doMultipliers ? Form("nu-tki-ptmuconstraint-%s-comps-pt-multiplier_ratio-2track_%s_%s.C",sample.c_str(),varx.c_str(),vary.c_str()) : Form("nu-tki-ptmuconstraint-%s-comps-pt_ratio-2track_%s_%s.C",sample.c_str(),varx.c_str(),vary.c_str()));
    
  }
  else{  
    gc1->Print(doMultipliers ? Form("nu-tki-ptmuconstraint_%s-comps-pt-multiplier_2track_%s_%s.eps",sample.c_str(),varx.c_str(),vary.c_str()) : Form("nu-tki-ptmuconstraint-%s-comps-pt_2track_%s_%s.eps",sample.c_str(),varx.c_str(),vary.c_str()));
    gc1->Print(doMultipliers ? Form("nu-tki-ptmuconstraint-%s-comps-pt-multiplier_2track_%s_%s.png",sample.c_str(),varx.c_str(),vary.c_str()) : Form("nu-tki-ptmuconstraint-%s-comps-pt_2track_%s_%s.png",sample.c_str(),varx.c_str(),vary.c_str()));
    gc1->Print(doMultipliers ? Form("nu-tki-ptmuconstraint-%s-comps-pt-multiplier_2track_%s_%s.C",sample.c_str(),varx.c_str(),vary.c_str()) : Form("nu-tki-ptmuconstraint-%s-comps-pt_2track_%s_%s.C",sample.c_str(),varx.c_str(),vary.c_str()));
  }
  TCanvas *cleg = new TCanvas("legend","legend",10,10,700,700);
  leg->Draw("");
  cleg->Print("Legend_bkg_samples_3D.png");
}

int main(int argc, char* argv[])
{
  
  string s_doRatio = argv[2];
  bool doRatio = s_doRatio=="1" ? true:false;
  cout << "Running with " << s_doRatio << "\t" << argv[3] << "\t" << argv[4] << endl;
  if(doRatio){
    makePlots(false,argv[1],"MichelSideBand",doRatio,argv[3],argv[4]);
    makePlots(false,argv[1],"BlobSideBand",doRatio,argv[3],argv[4]);
    makePlots(false,argv[1],"MicBlobSideBand",doRatio,argv[3],argv[4]);
  }
  else{
    makePlots(true,argv[1],"MichelSideBand",doRatio,argv[3],argv[4]);
    makePlots(true,argv[1],"BlobSideBand",doRatio,argv[3],argv[4]);
    makePlots(true,argv[1],"MicBlobSideBand",doRatio,argv[3],argv[4]);
  }

  return 0;
}
