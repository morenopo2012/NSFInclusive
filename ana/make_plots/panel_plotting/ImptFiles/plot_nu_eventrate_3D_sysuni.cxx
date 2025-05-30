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
//#include "Cintex/Cintex.h"

#include "myPlotStyle.h"
#include <algorithm>
#include "plot.h"

using namespace PlotUtils;

void makePlots(bool doMultipliers,  string location, bool doRatio, bool doZoom, string error_band)
{

  //ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);


  vector<double> recoil3Dbins;
  vector<double> pt3Dbins;
  vector<double> pz3Dbins;

  pt3Dbins.push_back(0.0);
  //  pt3Dbins.push_back(0.075);//added ME
  pt3Dbins.push_back(0.15);
  pt3Dbins.push_back(0.25);//added ME
  pt3Dbins.push_back(0.4);
  pt3Dbins.push_back(0.7);//added ME
  pt3Dbins.push_back(1.0);
  //  pt3Dbins.push_back(1.5);//added ME
  pt3Dbins.push_back(2.5);

  pz3Dbins.push_back(1.5);
  pz3Dbins.push_back(3.5);//added ME
  pz3Dbins.push_back(8.0);
  pz3Dbins.push_back(20.0);

  for(int i=0;i<10;i++)recoil3Dbins.push_back(i*0.04);
  for(int i=0;i<4;i++)recoil3Dbins.push_back(i*0.200+0.400);
  recoil3Dbins.push_back(1.9999);


  std::vector<std::vector<double> > full3D;
  full3D.push_back(recoil3Dbins);
  full3D.push_back(pt3Dbins);
  full3D.push_back(pz3Dbins);

  
  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);




  //three files 1-track, 2+track, N-track
  //CV
  TFile f1(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_CombinedPlaylists.root",location.c_str()));//1track
  TFile f2(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_CombinedPlaylists.root",location.c_str()));//2track
  TFile f3(Form("%s_CV/MuonEventSelection_MakeFlux-1_Multiplicity-0_Sample-Signal_CombinedPlaylists.root",location.c_str()));//Ntrack
  //Constraint File
  TFile f5(Form("%s_CV/SideBand3D_CombinedPlaylists.root",location.c_str()));

  //need pzmuptmu
  MnvH2D* mcMnv=(MnvH2D*)f3.Get("h_pzptrec_mc");//Get from N track
  MnvH2D* dataMnv = (MnvH2D*)f3.Get("h_pzptrec_data");//Get from N track
  
  MnvH2D* mcMnv_sys_uni0 = new MnvH2D(*mcMnv->GetVertErrorBand(error_band)->GetHist(0));
  MnvH2D* mcMnv_sys_uni1 = new MnvH2D(*mcMnv->GetVertErrorBand(error_band)->GetHist(1));



  std::vector<TH2D*> dataStat = my3d->Get2DHistos(dataMnv,false);
  std::vector<TH2D*> data = my3d->Get2DHistos(dataMnv,true);
  std::vector<TH2D*> mc = my3d->Get2DHistos(mcMnv,false);
  std::vector<TH2D*> mc_uni0 = my3d->Get2DHistos(mcMnv_sys_uni0,false);
  std::vector<TH2D*> mc_uni1 = my3d->Get2DHistos(mcMnv_sys_uni1,false);



  for(int i=1;i<pz3Dbins.size();i++){

    double pzwidth = pz3Dbins[i]-pz3Dbins[i-1];

    dataStat[i]->Scale(1e-5/pzwidth,"width");
    data[i]->Scale(1e-5/pzwidth,"width");
    mc[i]->Scale(1e-5/pzwidth,"width");
    mc_uni0[i]->Scale(1e-5/pzwidth,"width");
    mc_uni1[i]->Scale(1e-5/pzwidth,"width");

    if(doRatio){

    dataStat[i]->Divide(mc[i]);
    data[i]->Divide(mc[i]);
    mc_uni0[i]->Divide(mc[i]);
    mc_uni1[i]->Divide(mc[i]);
    mc[i]->Divide(mc[i]);

    }
    

    // These line and marker styles will be propagated to the 1D plots
    vector<int> mycolors = getColors(2);
    vector<int> mycolors2 = getColors(1);
    mc[i]->SetLineColor(kRed);
    mc[i]->SetLineWidth(2);
    
    mc_uni0[i]->SetLineColor(kViolet);
    mc_uni1[i]->SetLineColor(kBlue);

    
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
    histAndOpts.push_back(std::make_pair(mc[i],       "hist"));
    histAndOpts.push_back(std::make_pair(mc_uni0[i],       "hist"));
    histAndOpts.push_back(std::make_pair(mc_uni1[i],       "hist"));
    histAndOpts.push_back(std::make_pair(dataStat[i], "graph e"));
    histAndOpts.push_back(std::make_pair(data[i],     "graph ep"));
    
    
    
    // ----------------------------------------------------------------------------------
    //
    // First make pt in bins of pz
    
    // Values to multiply each bin by to get them on a similar range
    double multipliers[] ={1,1,1,1,1,1,2,2,2,2,2,5,10,20,1};
    double multipliers2[]={0.5,0.5,0.5,0.5,0.5,0.5,50,1,1,1,1,2,5,25,1};
    double multipliers3[]={10,10,10,10,10,10,10,20,20,20,20,20,100,300};
    double multipliers4[]={10,10,10,50,100,1,1,1,1,1,1,1,1,1,1};
    GridCanvas* gc=NULL;
    if(i==1) gc = plotYAxis1D(histAndOpts,"P_{t} (GeV)","Vis. E (GeV)", doMultipliers ? multipliers : NULL);
    if(i==2) gc = plotYAxis1D(histAndOpts,"P_{t} (GeV)","Vis. E (GeV)", doMultipliers ? multipliers2 : NULL);
    if(i==3) gc = plotYAxis1D(histAndOpts,"P_{t} (GeV)","Vis. E (GeV)", doMultipliers ? multipliers3 : NULL);
    if(i==4) gc = plotYAxis1D(histAndOpts,"P_{t} (GeV)","Vis. E (GeV)", doMultipliers ? multipliers4 : NULL);
    if(i==5) gc = plotYAxis1D(histAndOpts,"P_{t} (GeV)","Vis. E (GeV)", doMultipliers ? multipliers2 : NULL);
    // Set the y range manually. Can also use gc->Remax() to guess automatically
    if(doRatio) gc->SetYLimits(0.5,1.49);
    else gc->SetYLimits(0, 5);
    if(doZoom)gc->SetXLimits(0,0.4);
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
    leg->AddEntry(mc_uni0[i], "MINERvA Tune_sys_uni0", "l");
    leg->AddEntry(mc_uni1[i], "MINERvA Tune_sys_uni1", "l");

    TLatex mytex;
    mytex.SetTextSize(0.05);
    string mystring =     Form("%.2f < P_{||} [GeV] < %.2f",pz3Dbins[i-1],pz3Dbins[i]);
    mytex.DrawLatex(0.35,0.96,mystring.c_str());    
    leg->Draw("SAME");
    gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-model-pt-bin-%d.eps",i));
    gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-model-pt-bin-%d.png",i));
    gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-model-pt-bin-%d.C",i));
    
    // ----------------------------------------------------------------------------------
    //
    // Now make pz in bins of pt. It's all the same
    
    // Values to multiply each bin by to get them on a similar range
    double  multipliers5[]={4,1,1,1,4,200,75,1,1,1,1,1,1,1};
    double  multipliers6[]={2,0.7,0.5,0.5,1,10,5,500,1,1,1,1,1,1};
    double  multipliers7[]={25,10,10,10,20,100,5,100,1,1,1,1,1,1};
    double  multipliers8[]={100,50,25,25,25,50,100,500,1,1,1,1,1,1};
			    
    
    // plotXAxis1D fiddles the x axis values to squash up the tail so it
    // doesn't take up all the horizontal space.
    GridCanvas* gc2=NULL;
    if(i==1) gc2 = plotXAxis1D(histAndOpts,"Visible Energy (GeV)","P_{t} (GeV)", doMultipliers ? multipliers5 : NULL);
    if(i==2) gc2 = plotXAxis1D(histAndOpts,"Visible Energy (GeV)","P_{t} (GeV)", doMultipliers ? multipliers6 : NULL);
    if(i==3) gc2 = plotXAxis1D(histAndOpts,"Visible Energy (GeV)","P_{t} (GeV)", doMultipliers ? multipliers7 : NULL);
    if(i==4) gc2 = plotXAxis1D(histAndOpts,"Visible Energy (GeV)","P_{t} (GeV)", doMultipliers ? multipliers8 : NULL);
    if(i==5) gc2 = plotXAxis1D(histAndOpts,"Visible Energy (GeV)","P_{t} (GeV)", doMultipliers ? multipliers4 : NULL);

    mytex.DrawLatex(0.35,0.96,mystring.c_str());
    if(doRatio) gc2->SetYLimits(0.5,1.49);
    else gc2->SetYLimits(0, 5);
    if(doZoom)gc2->SetXLimits(0,0.4);
    gc2->SetYTitle("Event Rate(x10^{-5}) per GeV^{3}");
    gc2->Modified();
    gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-model-pz-bin-%d.eps",i));
    gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-model-pz-bin-%d.png",i));
    gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-model-pz-bin-%d.C",i));
  }
}

int main(int argc, char* argv[])
{

  string location = argv[1];
  string s_zoom = argv[2];
  string s_err = argv[3];
  bool doZoom = s_zoom=="1"? true:false;
  //multipliers
  //makePlots(true,location,false,doZoom,s_err);
  //standard
  makePlots(false,location,true,doZoom,s_err);

  return 0;
}
