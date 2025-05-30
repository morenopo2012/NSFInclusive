//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/HyperDimLinearizer.h"//THIS HAS TO CHANGE TO BE INCLUDED IN THE MAKE FILE EVENTUALLY.
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



void makePlots(bool doMultipliers, string location,string sample, bool doRatio)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  vector<double> recoil3Dbins;
  vector<double> pt3Dbins;
  vector<double> pz3Dbins;

  pt3Dbins.push_back(0.0);
  pt3Dbins.push_back(0.075);//added
  pt3Dbins.push_back(0.15);
  pt3Dbins.push_back(0.25);
  pt3Dbins.push_back(0.325);//added
  pt3Dbins.push_back(0.4);
  pt3Dbins.push_back(0.475);//added
  pt3Dbins.push_back(0.55);
  pt3Dbins.push_back(0.7);
  pt3Dbins.push_back(0.85);
  pt3Dbins.push_back(1.0);
  //  pt3Dbins.push_back(1.25);
  //  pt3Dbins.push_back(1.5);
  pt3Dbins.push_back(2.5);

  

  pz3Dbins.push_back(1.5);
  pz3Dbins.push_back(3.5);//added ME
  //These were added to fix unfolding
  pz3Dbins.push_back(4.5);//fix unf
  pz3Dbins.push_back(7.0);//fix unf
  pz3Dbins.push_back(8.0);
  //These are added to fix unfolding
  pz3Dbins.push_back(10.0);//fix unf
  pz3Dbins.push_back(20.0);

  recoil3Dbins.push_back(0.0);
  recoil3Dbins.push_back(20.0);
  //for(int i=0;i<10;i++)recoil3Dbins.push_back(i*40); 
  //for(int i=0;i<9;i++)recoil3Dbins.push_back(i*40+40);//40,80,120,160,200,240,280,320,360
  for(int i=0;i<4;i++)recoil3Dbins.push_back(i*40+40);//40,80,120,160
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*80+240);//240,320
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*200+400);
  recoil3Dbins.push_back(799.0);
  for(int i=0;i<recoil3Dbins.size();i++) recoil3Dbins[i]=recoil3Dbins[i]/1000.;


  std::vector<std::vector<double> > full3D;
  full3D.push_back(recoil3Dbins);
  full3D.push_back(pt3Dbins);
  full3D.push_back(pz3Dbins);

  
  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);



  TFile f1(Form("%s/MuonEventSelection_MakeFlux-1_Multiplicity-0_Sample-%s_CombinedPlaylists.root",location.c_str(),sample.c_str()));//sample 1 trk

  string scaleform = "";
  if(sample=="MichelSideBand") scaleform = "_michel";
  if(sample=="BlobSideBand") scaleform = "_blobs";
  if(sample=="MicBlobSideBand") scaleform = "_micblob";
  
  MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_recoil_ptmu_data");
  MnvH2D* mcMnv=(MnvH2D*)f1.Get("h_recoil_ptmu_mc");
  MnvH2D* mcsigMnv=(MnvH2D*)f1.Get("h_recoil_ptmu_qelike");
  MnvH2D* mcbkgMnv_scp=(MnvH2D*)f1.Get("h_recoil_ptmu_qelikenot_singlechargedpion");
  MnvH2D* mcbkgMnv_snp=(MnvH2D*)f1.Get("h_recoil_ptmu_qelikenot_singleneutralpion");
  MnvH2D* mcbkgMnv_mp=(MnvH2D*)f1.Get("h_recoil_ptmu_qelikenot_multipion");
  MnvH2D* mcbkgMnv_nop=(MnvH2D*)f1.Get("h_recoil_ptmu_qelikenot_not_scp_snp_mp");


  TH2* dataresults = new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mcresults = new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mcresults_qelike = new TH2D(mcsigMnv->GetCVHistoWithStatError());
  TH2* mcresults_scp = new TH2D(mcbkgMnv_scp->GetCVHistoWithStatError());
  TH2* mcresults_snp = new TH2D(mcbkgMnv_snp->GetCVHistoWithStatError());
  TH2* mcresults_mp = new TH2D(mcbkgMnv_mp->GetCVHistoWithStatError());
  TH2* mcresults_nop = new TH2D(mcbkgMnv_nop->GetCVHistoWithStatError());
  
  
  vector<int> mycolors = getColors(2);
  mcresults->SetLineColor(kRed);
  mcresults->SetLineWidth(2);
  
  dataresults->SetMarkerStyle(kFullCircle);
  dataresults->SetMarkerSize(0.7);
  dataresults->SetLineColor(kBlack);
  dataresults->SetLineWidth(2);

  mcresults_qelike->SetLineColor(kViolet);
  mcresults_scp->SetLineColor(kBlue);
  mcresults_snp->SetLineColor(kGreen);
  mcresults_mp->SetLineColor(kRed+3);
  mcresults_nop->SetLineColor(kBlack);
  /*
  dataresults->Scale(1,"width");
  mcresults_qelike->Scale(1,"width");
  mcresults_scp->Scale(1,"width");
  mcresults_snp->Scale(1,"width");
  mcresults_mp->Scale(1,"width");
  mcresults_nop->Scale(1,"width");
  mcresults->Scale(1,"width");
  */

  if(doRatio){
    dataresults->Divide(mcresults);
    mcresults_qelike->Divide(mcresults);
    mcresults_scp->Divide(mcresults);
    mcresults_snp->Divide(mcresults);
    mcresults_mp->Divide(mcresults);
    mcresults_nop->Divide(mcresults);
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
  histAndOpts.push_back(std::make_pair(mcresults,       "graph0 lx"));
  histAndOpts.push_back(std::make_pair(mcresults_qelike,       "graph0 lx"));
  histAndOpts.push_back(std::make_pair(mcresults_scp,       "graph0 lx"));
  histAndOpts.push_back(std::make_pair(mcresults_snp,       "graph0 lx"));
  histAndOpts.push_back(std::make_pair(mcresults_mp,       "graph0 lx"));
  histAndOpts.push_back(std::make_pair(mcresults_nop,       "graph0 lx"));
  
  cout << "Mults" << endl;
  

  vector<double> multi_1 = GetScales(histAndOpts, true);
  
    
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.79, 0.1, 1.0, 0.4);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(dataresults, "MINERvA data", "lpe");
  leg->AddEntry(mcresults, "Minerva Tune v1", "l");
  leg->AddEntry(mcresults_qelike, "QElike", "l");
  leg->AddEntry(mcresults_scp, "Charged Pion", "l");
  leg->AddEntry(mcresults_snp, "Neutral Pion", "l");
  leg->AddEntry(mcresults_mp, "Multi-Pion", "l");
  leg->AddEntry(mcresults_nop, "No Pion", "l");
  
  
  TLatex mytex;
  mytex.SetTextSize(0.05);
  //  string mystring =     Form("%.2f < P_{||} (GeV/c) < %.2f",pz3Dbins[i-1],pz3Dbins);
  
  
  GridCanvas* gc1=NULL;
  gc1=plotXAxis1D_ReducedXRange(histAndOpts, "Visible Energy (GeV)","P_{t}", 0,750,doMultipliers ? &multi_1[0] : NULL);
  
  // Set the y range manually. Can also use gc1->Remax() to guess automatically
  if(doRatio)gc1->SetYLimits(0,1.99);
  else gc1->SetYLimits(0, 100000);
  //gc1->SetYLimits(0, 1e-39);
  if(doRatio)gc1->SetYTitle("Ratio to Minerva Tune v1");
  else gc1->SetYTitle("Events per GeV^{3}/c^{2}");
  leg->Draw("SAME");
  gc1->Modified();
  
  if(doRatio){
    gc1->Print(doMultipliers ? Form("nu-3d-bkg-2dinputs_%s-comps-pt-multiplier_-ratio-alltrack.eps",sample.c_str()) : Form("nu-3d-bkg-2dinputs-%s-comps-pt_-ratio-alltrack.eps",sample.c_str()));
    gc1->Print(doMultipliers ? Form("nu-3d-bkg-2dinputs-%s-comps-pt-multiplier_-ratio-alltrack.png",sample.c_str()) : Form("nu-3d-bkg-2dinputs-%s-comps-pt_-ratio-alltrack.png",sample.c_str()));
    gc1->Print(doMultipliers ? Form("nu-3d-bkg-2dinputs-%s-comps-pt-multiplier_-ratio-alltrack.C",sample.c_str()) : Form("nu-3d-bkg-2dinputs-%s-comps-pt_-ratio-alltrack.C",sample.c_str()));
  }
  else{
    gc1->Print(doMultipliers ? Form("nu-3d-bkg-2dinputs_%s-comps-pt-multiplier_-alltrack.eps",sample.c_str()) : Form("nu-3d-bkg-2dinputs-%s-comps-pt_-alltrack.eps",sample.c_str()));
    gc1->Print(doMultipliers ? Form("nu-3d-bkg-2dinputs-%s-comps-pt-multiplier_-alltrack.png",sample.c_str()) : Form("nu-3d-bkg-2dinputs-%s-comps-pt_-alltrack.png",sample.c_str()));
    gc1->Print(doMultipliers ? Form("nu-3d-bkg-2dinputs-%s-comps-pt-multiplier_-alltrack.C",sample.c_str()) : Form("nu-3d-bkg-2dinputs-%s-comps-pt_-alltrack.C",sample.c_str()));
  }
  TCanvas *cleg = new TCanvas("legend","legend",10,10,700,700);
  leg->Draw("");
  cleg->Print("Legend_bkg_samples_3D.png");
}

int main(int argc, char* argv[])
{
  
  string s_doRatio = argv[2];
  bool doRatio = s_doRatio=="1" ? true:false;
  if(doRatio){
    makePlots(false,argv[1],"MichelSideBand",doRatio);
    makePlots(false,argv[1],"BlobSideBand",doRatio);
    makePlots(false,argv[1],"MicBlobSideBand",doRatio);
    makePlots(false,argv[1],"Signal",doRatio);
  }
  else{
    makePlots(false,argv[1],"MichelSideBand",doRatio);
    makePlots(false,argv[1],"BlobSideBand",doRatio);
    makePlots(false,argv[1],"MicBlobSideBand",doRatio);
    makePlots(false,argv[1],"Signal",doRatio);
  }

  return 0;
}
