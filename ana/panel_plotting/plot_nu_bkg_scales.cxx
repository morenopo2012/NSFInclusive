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



void makePlots(bool doMultipliers, string location, string location2)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  MnvPlotter *plotter = new MnvPlotter();
  TFile f1(Form("%s",location.c_str()));//sample 1 trk
  TFile f2(Form("%s",location2.c_str()));//sample 1 trk
  f1.ls();
  

  MnvH2D* scalefactorNotQE = (MnvH2D*)f1.Get("Signal_qelikenot_scaled");
  MnvH2D* scalefactorQE = (MnvH2D*)f1.Get("qescale");
  MnvH2D* scalefactorSCP = (MnvH2D*)f1.Get("scpscale");
  MnvH2D* scalefactorSNP = (MnvH2D*)f1.Get("snpscale");
  MnvH2D* scalefactorMP = (MnvH2D*)f1.Get("mpscale");

  TH2D* sfNotQE = new TH2D(scalefactorNotQE->GetCVHistoWithStatError());
  TH2D* sfQE = new TH2D(scalefactorQE->GetCVHistoWithStatError());
  TH2D* sfSCP = new TH2D(scalefactorSCP->GetCVHistoWithStatError());
  TH2D* sfSNP = new TH2D(scalefactorSNP->GetCVHistoWithStatError());
  TH2D* sfMP = new TH2D(scalefactorMP->GetCVHistoWithStatError());


  MnvH2D* scalefactorNotQE2 = NULL;
  TH2D* sfNotQE2 = NULL;
  TH2D* sfNotQERatio = NULL;
  bool doSFComp = location2!="";
  if(doSFComp){
    scalefactorNotQE2 = (MnvH2D*)f2.Get("Signal_qelikenot_scaled");
    sfNotQE2  = new TH2D(scalefactorNotQE2->GetCVHistoWithStatError());
    sfNotQERatio = (TH2D*)sfNotQE2->Clone("Ratio");
    sfNotQERatio->Divide(sfNotQE);
    
  }


  
  vector<int> mycolors = getColors(2);

  sfNotQE->SetLineColor(kBlack);
  sfNotQE->SetLineWidth(3);
  if(doSFComp){
    sfNotQE2->SetLineStyle(2);
    sfNotQE2->SetLineWidth(2);
    sfNotQERatio->SetLineColor(kRed);
  }
  sfQE->SetLineColor(kRed);
  sfSCP->SetLineColor(mycolors[3]);
  sfSNP->SetLineColor(mycolors[4]);
  sfMP->SetLineColor(mycolors[5]);
   

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
  if(!doSFComp){
    histAndOpts.push_back(std::make_pair(sfNotQE, "histl"));
    histAndOpts.push_back(std::make_pair(sfQE, "histl"));
    histAndOpts.push_back(std::make_pair(sfSCP, "histl"));
    histAndOpts.push_back(std::make_pair(sfSNP, "histl"));
    histAndOpts.push_back(std::make_pair(sfMP, "histl"));
  }
  else{
    histAndOpts.push_back(std::make_pair(sfNotQE, "histl"));
    histAndOpts.push_back(std::make_pair(sfNotQE2, "histl"));
    histAndOpts.push_back(std::make_pair(sfNotQERatio, "histl"));
  }
  vector<double> multi_1 = GetScales(histAndOpts, true);
  
    
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.1, 0.91, 1, 1);
  leg->SetNColumns(3);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  if(!doSFComp){
    leg->AddEntry(sfNotQE, "Combined Bkg Constraint","L");
    leg->AddEntry(sfQE, "QE-Like Scale","L");
    leg->AddEntry(sfSCP, "Single Charged Pion Scale","L");
    leg->AddEntry(sfSNP, "Single Neutral Pion Scale", "L");
    leg->AddEntry(sfMP, "Multi-pion Scale","L");
  }
  else{
    leg->AddEntry(sfNotQE, "Combined Bkg Constraint","L");
    leg->AddEntry(sfNotQE2, "Combined Bkg Constraint 2","L");
    leg->AddEntry(sfNotQERatio, "Ratio","L");
  }
  
  
  
  GridCanvas* gc1=NULL;
  gc1=plotXAxis1D_ReducedXRange(histAndOpts, "Visible Energy (GeV)","P_{t}", 0,1000,doMultipliers ? &multi_1[0] : NULL);
  //  gc1->SetLogx(true);
  // Set the y range manually. Can also use gc1->Remax() to guess automatically
  gc1->SetYLimits(-1, 3.99);
  //gc1->SetYLimits(0, 1e-39);
  gc1->SetYTitle("Sub-sample scale factor");
  leg->Draw("SAME");
  gc1->Modified();
  
  if(!doSFComp){
    gc1->Print(doMultipliers ? "nu-3d-bkg-2dscales-comps-pt-multiplier_alltrack.eps" : "nu-3d-bkg-2dscales-comps-pt_altrack.eps");
    gc1->Print(doMultipliers ? "nu-3d-bkg-2dscales-comps-pt-multiplier_alltrack.png" : "nu-3d-bkg-2dscales-comps-pt_alltrack.png");
    gc1->Print(doMultipliers ? "nu-3d-bkg-2dscales-comps-pt-multiplier_alltrack.C" : "nu-3d-bkg-2dscales-comps-pt_alltrack.C");
  }
  else{
    gc1->Print(doMultipliers ? "nu-3d-bkg-2dscales_compTotalBkg-comps-pt-multiplier_alltrack.eps" : "nu-3d-bkg-2dscales-comps-pt_altrack.eps");
    gc1->Print(doMultipliers ? "nu-3d-bkg-2dscales_compTotalBkg-comps-pt-multiplier_alltrack.png" : "nu-3d-bkg-2dscales-comps-pt_alltrack.png");
    gc1->Print(doMultipliers ? "nu-3d-bkg-2dscales_compTotalBkg-comps-pt-multiplier_alltrack.C" : "nu-3d-bkg-2dscales-comps-pt_alltrack.C");
  }

  GridCanvas* gc2=NULL;
  gc2=plotXAxis1D_ReducedXRange(histAndOpts, "Visible Energy (GeV)","P_{t}", 0,40000,doMultipliers ? &multi_1[0] : NULL);
  gc2->SetLogx(true);
  // Set the y range manually. Can also use gc2->Remax() to guess automatically
  gc2->SetYLimits(-4, 3.99);
  //gc2->SetYLimits(0, 1e-39);
  gc2->SetYTitle("Sub-sample scale factor");
  leg->Draw("SAME");
  gc2->Modified();
  if(!doSFComp){
    gc2->Print(doMultipliers ? "nu-3d-bkg-2dscales-comps-pt-multiplier_allrange_alltrack.eps" : "nu-3d-bkg-2dscales-comps-pt_altrack.eps");
    gc2->Print(doMultipliers ? "nu-3d-bkg-2dscales-comps-pt-multiplier_allrange_alltrack.png" : "nu-3d-bkg-2dscales-comps-pt_allrange_alltrack.png");
    gc2->Print(doMultipliers ? "nu-3d-bkg-2dscales-comps-pt-multiplier_allrange_alltrack.C" : "nu-3d-bkg-2dscales-comps-pt_allrange_alltrack.C");
  }
  else{
    gc2->Print(doMultipliers ? "nu-3d-bkg-2dscales_compTotalBkg-comps-pt-multiplier_allrange_alltrack.eps" : "nu-3d-bkg-2dscales-comps-pt_altrack.eps");
    gc2->Print(doMultipliers ? "nu-3d-bkg-2dscales_compTotalBkg-comps-pt-multiplier_allrange_alltrack.png" : "nu-3d-bkg-2dscales-comps-pt_allrange_alltrack.png");
    gc2->Print(doMultipliers ? "nu-3d-bkg-2dscales_compTotalBkg-comps-pt-multiplier_allrange_alltrack.C" : "nu-3d-bkg-2dscales-comps-pt_allrange_alltrack.C");
  }

  GridCanvas* gc3=NULL;
  gc3=plotXAxis1D_ReducedXRange(histAndOpts, "Visible Energy (GeV)","P_{t}", 0,40000,doMultipliers ? &multi_1[0] : NULL);
  gc3->SetLogx(true);
  // Set the y range manually. Can also use gc3->Remax() to guess automatically
  gc3->SetYLimits(0.75, 1.25);
  //gc3->SetYLimits(0, 1e-39);
  gc3->SetYTitle("Sub-sample scale factor");
  leg->Draw("SAME");
  gc3->Modified();
  if(!doSFComp){
    gc3->Print(doMultipliers ? "nu-3d-bkg-2dscales-comps-pt-multiplier_zoomy_allrange_alltrack.eps" : "nu-3d-bkg-2dscales-comps-pt_altrack.eps");
    gc3->Print(doMultipliers ? "nu-3d-bkg-2dscales-comps-pt-multiplier_zoomy_allrange_alltrack.png" : "nu-3d-bkg-2dscales-comps-pt_zoomy_allrange_alltrack.png");
    gc3->Print(doMultipliers ? "nu-3d-bkg-2dscales-comps-pt-multiplier_zoomy_allrange_alltrack.C" : "nu-3d-bkg-2dscales-comps-pt_zoomy_allrange_alltrack.C");
  }
  else{
    gc3->Print(doMultipliers ? "nu-3d-bkg-2dscales_compTotalBkg-comps-pt-multiplier_zoomy_allrange_alltrack.eps" : "nu-3d-bkg-2dscales-comps-pt_altrack.eps");
    gc3->Print(doMultipliers ? "nu-3d-bkg-2dscales_compTotalBkg-comps-pt-multiplier_zoomy_allrange_alltrack.png" : "nu-3d-bkg-2dscales-comps-pt_zoomy_allrange_alltrack.png");
    gc3->Print(doMultipliers ? "nu-3d-bkg-2dscales_compTotalBkg-comps-pt-multiplier_zoomy_allrange_alltrack.C" : "nu-3d-bkg-2dscales-comps-pt_zoomy_allrange_alltrack.C");
  }
  
}

int main(int argc, char* argv[])
{
  
  makePlots(false,argv[1],"");
  if(argc==3) makePlots(false,argv[1],argv[2]);

  return 0;
}
