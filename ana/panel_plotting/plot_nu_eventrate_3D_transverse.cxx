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
//#include "localColor.h"
#include "PlotUtils/MnvColors.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;

void makePlots(bool doMultipliers, string location,  bool doRatio, string varx, string varz)
{
  if(doRatio) doMultipliers=false;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  map<string,vector<double> > binning_by_var;


  vector<double> recoil3Dbins;
  recoil3Dbins.push_back(0.0);
  recoil3Dbins.push_back(0.020);
  for(int i=0;i<4;i++)recoil3Dbins.push_back(i*0.040+0.040);//40,80,120,160
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*0.080+0.240);//240,320
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*0.200+0.400);
  recoil3Dbins.push_back(0.799);

  vector<double> pt3Dbins;
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
  pt3Dbins.push_back(1.25);
  pt3Dbins.push_back(1.5);
  pt3Dbins.push_back(2.5);
  vector<double> pz3Dbins;
  pz3Dbins.push_back(1.5);
  pz3Dbins.push_back(3.5);//added ME
  pz3Dbins.push_back(4.5);
  pz3Dbins.push_back(7.0);
  pz3Dbins.push_back(8.0);
  pz3Dbins.push_back(10.0);
  pz3Dbins.push_back(20.0);
  vector<double> dalphatbins;
  for(int i=0; i<8;i++) dalphatbins.push_back(i*15);
  for(int i=0; i<7;i++) dalphatbins.push_back(i*10+120);
  vector<double> dphitbins;
  dphitbins.push_back(0);
  dphitbins.push_back(5);
  dphitbins.push_back(10);
  dphitbins.push_back(15);
  dphitbins.push_back(20);
  dphitbins.push_back(25);
  dphitbins.push_back(30);
  dphitbins.push_back(40);
  dphitbins.push_back(50);
  dphitbins.push_back(60);
  dphitbins.push_back(70);
  dphitbins.push_back(80);
  dphitbins.push_back(90);
  dphitbins.push_back(100);
  dphitbins.push_back(110);
  dphitbins.push_back(120);
  dphitbins.push_back(130);
  dphitbins.push_back(140);
  dphitbins.push_back(150);
  dphitbins.push_back(160);
  dphitbins.push_back(170);
  dphitbins.push_back(180);
  vector<double> dptbins;
  for( int i = 0; i < 20; i++ ){
    dptbins.push_back(i*0.05);
  }
  for( int i = 0; i<5; i++){
    dptbins.push_back(i*0.2+1);
  }
  dptbins.push_back(2);
  dptbins.push_back(2.5);
  vector<double> dptxbins;
  dptxbins.push_back(-2.0);
  //  dptxbins.push_back(-1.0);
  dptxbins.push_back(-0.5);
  dptxbins.push_back(-0.4);
  dptxbins.push_back(-0.3);
  for(int i=0;i<5;i++) dptxbins.push_back(i*0.1-0.2);
  dptxbins.push_back(0.3);
  dptxbins.push_back(0.4);
  dptxbins.push_back(0.5);
  //  dptxbins.push_back(1.0);
  dptxbins.push_back(2.0);
  vector<double> dptybins;
  dptybins.push_back(-2.0);
  dptybins.push_back(-1.30);
  for(int i=0;i<15;i++) dptybins.push_back(i*0.1-1.2);
  dptybins.push_back(0.3);
  dptybins.push_back(0.4);
  dptybins.push_back(0.5);
  dptybins.push_back(1.0);
  vector<double> pnbins;
  for( int i = 0; i < 20; i++ ){
    pnbins.push_back(i*0.05);
  }
  for( int i = 0; i<5; i++){
    pnbins.push_back(i*0.2+1);
  }
  pnbins.push_back(2);
  pnbins.push_back(2.5);
  vector<double> pthetabins;

  vector<double> pkinbins;
  pkinbins.push_back(0);
  pkinbins.push_back(0.1);
  pkinbins.push_back(0.125);//500 MeV/c Jeffrey TKi cuttoff
  pkinbins.push_back(0.14);
  pkinbins.push_back(0.16);
  pkinbins.push_back(0.18);
  pkinbins.push_back(0.20);
  pkinbins.push_back(0.22);
  pkinbins.push_back(0.24);
  pkinbins.push_back(0.26);
  pkinbins.push_back(0.28);
  pkinbins.push_back(0.30);
  pkinbins.push_back(0.32);
  pkinbins.push_back(0.34);
  pkinbins.push_back(0.36);
  pkinbins.push_back(0.38);
  pkinbins.push_back(0.40);
  pkinbins.push_back(0.42);
  pkinbins.push_back(0.44);
  pkinbins.push_back(0.46);
  pkinbins.push_back(0.48);
  pkinbins.push_back(0.51);
  pkinbins.push_back(0.54);
  pkinbins.push_back(0.585);//My 1.2 GeV/c cutoff
  vector<double> dplbins;
  dplbins.push_back(-0.7);
  dplbins.push_back(0);
  dplbins.push_back(0.08);
  dplbins.push_back(0.12);
  dplbins.push_back(0.16);
  dplbins.push_back(0.19);
  dplbins.push_back(0.22);
  dplbins.push_back(0.25);
  dplbins.push_back(0.28);
  dplbins.push_back(0.31);
  dplbins.push_back(0.34);
  dplbins.push_back(0.36);
  dplbins.push_back(0.38);
  dplbins.push_back(0.41);
  dplbins.push_back(0.44);
  dplbins.push_back(0.47);
  dplbins.push_back(0.52);
  dplbins.push_back(0.6);
  vector<double> thetapnbins;
  thetapnbins.push_back(-1.7);
  thetapnbins.push_back(-0.68);
  thetapnbins.push_back(-0.16);
  thetapnbins.push_back(0.16);
  thetapnbins.push_back(0.4);
  thetapnbins.push_back(0.6);
  thetapnbins.push_back(0.72);
  thetapnbins.push_back(0.8);
  thetapnbins.push_back(0.88);
  thetapnbins.push_back(0.92);
  thetapnbins.push_back(0.96);
  thetapnbins.push_back(1.0);
  thetapnbins.push_back(1.04);
  thetapnbins.push_back(1.08);
  thetapnbins.push_back(1.5);

  vector<double> signedbins;
  signedbins.push_back(0);
  signedbins.push_back(1);
  signedbins.push_back(2);


  binning_by_var["muonPz"] = pz3Dbins;
  binning_by_var["dalphat"] = dalphatbins;
  binning_by_var["dphit"] = dphitbins;
  binning_by_var["pn"] = pnbins;
  binning_by_var["dpt"] = dptbins;
  binning_by_var["dptx"] = dptxbins;
  binning_by_var["dpty"] = dptybins;
  binning_by_var["signed"] = signedbins;
  binning_by_var["recoil"]= recoil3Dbins;
  binning_by_var["dpl"] = dplbins;
  binning_by_var["thetapn"]=thetapnbins;
  binning_by_var["tp"]=pkinbins;






  std::vector<std::vector<double> > full3D;
  full3D.push_back(binning_by_var[varx]);
  full3D.push_back(pt3Dbins);
  full3D.push_back(binning_by_var[varz]);

  
  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);





  TFile f3(Form("%s",location.c_str()));//Ntrack


  //need pzmuptmu
  MnvH2D* mcMnv=(MnvH2D*)f3.Get(Form("h_%s_ptmu_%s_mc",varx.c_str(),varz.c_str()));//Get from N track
  MnvH2D* dataMnv = (MnvH2D*)f3.Get(Form("h_%s_ptmu_%s_data",varx.c_str(),varz.c_str()));//Get from N track

  MnvH2D* mcMnv_qe = (MnvH2D*)f3.Get(Form("h_%s_ptmu_%s_qelike_qe",varx.c_str(),varz.c_str()));//Get from N track
  MnvH2D* mcMnv_res = (MnvH2D*)f3.Get(Form("h_%s_ptmu_%s_qelike_res",varx.c_str(),varz.c_str()));//Get from N track
  MnvH2D* mcMnv_dis = (MnvH2D*)f3.Get(Form("h_%s_ptmu_%s_qelike_dis",varx.c_str(),varz.c_str()));//Get from N track
  MnvH2D* mcMnv_2p2h = (MnvH2D*)f3.Get(Form("h_%s_ptmu_%s_qelike_2p2h",varx.c_str(),varz.c_str()));//Get from N track
  MnvH2D* mcMnv_bkg = (MnvH2D*)f3.Get(Form("h_%s_ptmu_%s_qelikenot",varx.c_str(),varz.c_str()));//Get from N track
  

  std::vector<TH2D*> dataStat = my3d->Get2DHistos(dataMnv,false);
  std::vector<TH2D*> data = my3d->Get2DHistos(dataMnv,true);
  std::vector<TH2D*> mc = my3d->Get2DHistos(mcMnv,false);
  std::vector<TH2D*> mc_qe = my3d->Get2DHistos(mcMnv_qe,false);
  std::vector<TH2D*> mc_res = my3d->Get2DHistos(mcMnv_res,false);
  std::vector<TH2D*> mc_dis = my3d->Get2DHistos(mcMnv_dis,false);
  std::vector<TH2D*> mc_2p2h = my3d->Get2DHistos(mcMnv_2p2h,false);
  std::vector<TH2D*> mc_bkg = my3d->Get2DHistos(mcMnv_bkg,false);



  for(int i=1;i<binning_by_var[varz].size();i++){

    double pzwidth = binning_by_var[varz][i]-binning_by_var[varz][i-1];

    dataStat[i]->Scale(1e-5/pzwidth,"width");
    data[i]->Scale(1e-5/pzwidth,"width");
    mc[i]->Scale(1e-5/pzwidth,"width");
    mc_qe[i]->Scale(1e-5/pzwidth,"width");
    mc_res[i]->Scale(1e-5/pzwidth,"width");
    mc_dis[i]->Scale(1e-5/pzwidth,"width");
    mc_2p2h[i]->Scale(1e-5/pzwidth,"width");
    mc_bkg[i]->Scale(1e-5/pzwidth,"width");


    if(doRatio){

    dataStat[i]->Divide(mc[i]);
    data[i]->Divide(mc[i]);
    mc_qe[i]->Divide(mc[i]);
    mc_res[i]->Divide(mc[i]);
    mc_dis[i]->Divide(mc[i]);
    mc_2p2h[i]->Divide(mc[i]);
    mc_bkg[i]->Divide(mc[i]);
    mc[i]->Divide(mc[i]);

    }
    

    // These line and marker styles will be propagated to the 1D plots
    vector<int> mycolors = MnvColors::GetColors(9);
    vector<int> mycolors2 = MnvColors::GetColors(6);
    mc[i]->SetLineColor(kRed);
    mc[i]->SetLineWidth(2);
    
    //need to add signal and bkg colors
    mc_qe[i]->SetLineColor(mycolors[2]);
    mc_res[i]->SetLineColor(mycolors[1]);
    mc_dis[i]->SetLineColor(mycolors[7]);
    mc_2p2h[i]->SetLineColor(mycolors[0]);
    mc_bkg[i]->SetLineColor(mycolors2[5]);
    //need to add 1track and 2 track
    
    // These line and marker styles will be propagated to the 1D plots
    data[i]->SetMarkerStyle(kFullCircle);
    data[i]->SetMarkerSize(0.7);
    data[i]->SetLineColor(kBlack);
    data[i]->SetLineWidth(2);
    
    dataStat[i]->SetMarkerStyle(1);
    dataStat[i]->SetLineColor(kBlack);
    dataStat[i]->SetLineWidth(2);
    
    
    
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
    histAndOpts.push_back(std::make_pair(dataStat[i], "histp"));
    histAndOpts.push_back(std::make_pair(mc[i],       "hist"));
    
    histAndOpts.push_back(std::make_pair(mc_qe[i],       "hist"));
    histAndOpts.push_back(std::make_pair(mc_res[i],       "hist"));
    histAndOpts.push_back(std::make_pair(mc_dis[i],       "hist"));
    histAndOpts.push_back(std::make_pair(mc_2p2h[i],       "hist"));
    histAndOpts.push_back(std::make_pair(mc_bkg[i],       "hist"));
    histAndOpts.push_back(std::make_pair(dataStat[i], "graph e"));
    histAndOpts.push_back(std::make_pair(data[i],     "graph ep"));
    
    
    
    // ----------------------------------------------------------------------------------
    //
    // First make pt in bins of pz
    
    // Values to multiply each bin by to get them on a similar range
    vector<double> multipliers = GetScales(histAndOpts, false, 5,0.75);
    GridCanvas* gc = plotYAxis1D(histAndOpts,"P_{t} (GeV)",varx, doMultipliers ? &multipliers[0] : NULL);
    // Set the y range manually. Can also use gc->Remax() to guess automatically
    if(doRatio) gc->SetYLimits(0,1.99);
    else gc->SetYLimits(0, 5);
    gc->SetYTitle("Event Rate(x10^{-5}) per GeV^{3}");
    gc->Modified();
    // Example of adding a legend. The co-ordinate system is NDC on the
    // entire canvas, ie (0,0) in the bottom left corner of the canvas
    // (not the individual pad), and (1,1) in the top right
    TLegend* leg=new TLegend(0.7, 0.1, 0.9, 0.3);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(data[i], "MINERvA data", "lpe");
    leg->AddEntry(mc[i], "MINERvA Tune", "l");
    leg->AddEntry(mc_qe[i],"QE","l");
    leg->AddEntry(mc_res[i],"Resonant","l");
    leg->AddEntry(mc_dis[i],"DIS","l");    
    leg->AddEntry(mc_2p2h[i],"2p2h","l");
    leg->AddEntry(mc_bkg[i],"Backgrounds","l");
        
    

    TLatex mytex;
    mytex.SetTextSize(0.05);
    string mystring =     Form("%.2f < %s < %.2f",binning_by_var[varz][i-1],varz.c_str(),binning_by_var[varz][i]);
    mytex.DrawLatex(0.35,0.96,mystring.c_str());    
    leg->Draw("SAME");
    gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-model-pt-bin-%d.eps",i));
    gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-model-pt-bin-%d.png",i));
    gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-model-pt-bin-%d.C",i));
  
    
    // ----------------------------------------------------------------------------------
    //
    // Now make pz in bins of pt. It's all the same
    
    // Values to multiply each bin by to get them on a similar range
    vector<double> multipliers2 = GetScales(histAndOpts, true, 5,0.75);
    
    // plotXAxis1D fiddles the x axis values to squash up the tail so it
    // doesn't take up all the horizontal space.
    GridCanvas* gc2 = plotXAxis1D_ReducedXRange(histAndOpts,varx,"P_{t} (GeV)", 0,2.0,doMultipliers ? &multipliers2[0] : NULL);

    mytex.DrawLatex(0.35,0.96,mystring.c_str());
    if(doRatio) gc2->SetYLimits(0,1.99);
    else gc2->SetYLimits(0, 5);
    gc2->SetYTitle("Event Rate(x10^{-5}) per GeV^{3}");
    leg->Draw("SAME");
    gc2->Modified();
    gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-model-pz-bin-%d.eps",i));
    gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-model-pz-bin-%d.png",i));
    gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-model-pz-bin-%d.C",i));
  }
}

int main(int argc, char* argv[])
{

  string location = argv[1];
  string s_ratio = argv[2];

  bool doRatio = s_ratio=="1" ? true:false;
  //multipliers
  makePlots(true,location,doRatio,argv[3],argv[4]);


  return 0;
}
