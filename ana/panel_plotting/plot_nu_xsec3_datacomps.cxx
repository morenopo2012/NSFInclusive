
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

using namespace PlotUtils;
void makePlots(string location1, string location2, string sample1, string sample2, string histname, string histnamemc)
{

  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  
  TFile f1(location1.c_str());
  TFile f2(location2.c_str());

  MnvH2D* data1Mnv=(MnvH2D*)f1.Get(histname.c_str());
  MnvH2D* data2Mnv=(MnvH2D*)f2.Get(histname.c_str());
  MnvH2D* mc1Mnv=(MnvH2D*)f1.Get(histnamemc.c_str());
  MnvH2D* mc2Mnv=(MnvH2D*)f2.Get(histnamemc.c_str());


  vector<double> recoil3Dbins;
  vector<double> pt3Dbins;
  vector<double> pz3Dbins;

  //Finer version
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
  //
  pz3Dbins.push_back(8.0);
  //These are added to fix unfolding
  pz3Dbins.push_back(10.0);//fix unf
  pz3Dbins.push_back(20.0);


  recoil3Dbins.push_back(0.0);
  recoil3Dbins.push_back(20.0);
  for(int i=0;i<4;i++)recoil3Dbins.push_back(i*40+40);//40,80,120,160
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*80+240);//240,320
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*200+400);
  recoil3Dbins.push_back(799.0);
  for(int i=0;i<recoil3Dbins.size();i++) recoil3Dbins[i]=recoil3Dbins[i]/1000.;


  /*
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
  //These were added to fix unfolding
  pz3Dbins.push_back(4.5);//fix unf
  pz3Dbins.push_back(7.0);//fix unf
  //
  pz3Dbins.push_back(8.0);
  //These are added to fix unfolding
  pz3Dbins.push_back(10.0);//fix unf
  pz3Dbins.push_back(20.0);

  for(int i=0;i<10;i++)recoil3Dbins.push_back(i*0.04);
  for(int i=0;i<3;i++)recoil3Dbins.push_back(i*0.2+0.4);
  */
  /*
  for(int i=0;i<4;i++)recoil3Dbins.push_back(i*0.2+0.4);
  recoil3Dbins.push_back(1.9999);
  */

  std::vector<std::vector<double> > full3D;
  full3D.push_back(recoil3Dbins);
  full3D.push_back(pt3Dbins);
  full3D.push_back(pz3Dbins);
  
  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);

  std::cout << "Starting up getting the projections" << std::endl;
  std::vector<TH2D*> data1results = my3d->Get2DHistos(data1Mnv,true);
  std::vector<TH2D*> data2results = my3d->Get2DHistos(data2Mnv,true);


  MnvH2D *data1Mnv_clone = (MnvH2D*)data1Mnv->Clone("data1clone");
  MnvH2D *data2Mnv_clone = (MnvH2D*)data2Mnv->Clone("data2clone");

  //kill stat
  for(int l=0;l<data1Mnv_clone->GetNbinsX()+2;l++){
    for(int k=0;k<data1Mnv_clone->GetNbinsY()+2;k++){
      data1Mnv_clone->SetBinError(l,k,0);
      data2Mnv_clone->SetBinError(l,k,0);
      mc1Mnv->SetBinError(l,k,0);
      mc2Mnv->SetBinError(l,k,0);
    }
  }

  //remove other error bands
  vector<string> verts = data1Mnv_clone->GetVertErrorBandNames();
  vector<string> lats = data1Mnv_clone->GetLatErrorBandNames();
  data1Mnv_clone->PopVertErrorBand("Flux");
  data2Mnv_clone->PopVertErrorBand("Flux");
  for(int i=0;i<lats.size();i++){
    data1Mnv_clone->PopLatErrorBand(lats[i]);
    data2Mnv_clone->PopLatErrorBand(lats[i]);
  }
  
  


  std::vector<TH2D*> data1results_remove = my3d->Get2DHistos(data1Mnv_clone,true);
  std::vector<TH2D*> data2results_remove = my3d->Get2DHistos(data2Mnv_clone,true);

  std::vector<TH2D*> mc1results = my3d->Get2DHistos(mc1Mnv,false);
  std::vector<TH2D*> mc2results = my3d->Get2DHistos(mc2Mnv,false);


  for(unsigned int i=1;i<pz3Dbins.size();i++){
    double pzwidth = pz3Dbins[i]-pz3Dbins[i-1];
    data1results[i]->Scale(1e39/pzwidth,"width");
    data2results[i]->Scale(1e39/pzwidth,"width");
    mc1results[i]->Scale(1e39/pzwidth,"width");
    mc2results[i]->Scale(1e39/pzwidth,"width");
  }


  vector<int> mycolors = getColors(2);
  vector<int> mycolors2 = getColors(1);
  for(unsigned int i=1;i<pz3Dbins.size();i++){//skip underflow and overflow pz bins
    cout << "DOING BIN " << i << endl;
   
    TH2D *errorhist_scaled = (TH2D*)data2results_remove[i]->Clone("errorhist");
    data2results_remove[i]->SaveAs(Form("Hist_%d.root",i));
    for(int l=0;l<errorhist_scaled->GetNbinsX()+2;l++){
      for(int k=0;k<errorhist_scaled->GetNbinsY()+2;k++){
	double val = errorhist_scaled->GetBinContent(l,k);
	double valerr = errorhist_scaled->GetBinError(l,k);
	cout <<i << "\t" <<l << "\t" << k << "\t" << val << "\t" << valerr << "\t" << valerr/val << endl;
	errorhist_scaled->SetBinContent(l,k,1);
	errorhist_scaled->SetBinError(l,k,valerr/val);
      }
    }
    errorhist_scaled->SetFillColor(2);
    //    errorhist_scaled->SetFillStyle(3004);

    data1results[i]->Divide(mc2results[i]);
    data2results[i]->Divide(mc2results[i]);
    mc1results[i]->Divide(mc2results[i]);
    mc2results[i]->Divide(mc2results[i]);
    

    TH2D* dataratio = (TH2D*)data1results[i]->Clone("Model");
    dataratio->Divide(data2results[i]);
    dataratio->SetMarkerStyle(8);
    dataratio->SetLineColor(1);
    dataratio->SetMarkerColor(1);

    data1results[i]->SetMarkerColor(mycolors[2]);
    mc1results[i]->SetLineColor(mycolors[2]);
    data2results[i]->SetMarkerColor(kBlack);
    data2results[i]->SetLineColor(1);
    mc2results[i]->SetLineColor(2);
    data1results[i]->SetMarkerStyle(8);
    data2results[i]->SetMarkerStyle(3);

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
    histAndOpts.push_back(std::make_pair(data1results[i], "histel"));
    histAndOpts.push_back(std::make_pair(data2results[i], "histpe"));
    histAndOpts.push_back(std::make_pair(mc1results[i], "hist"));
    histAndOpts.push_back(std::make_pair(mc2results[i], "hist"));

    std::vector<std::pair<TH2*, const char*> > histAndOpts2;
    histAndOpts2.push_back(std::make_pair(dataratio, "histl"));
    histAndOpts2.push_back(std::make_pair(errorhist_scaled, "pe3"));
    histAndOpts2.push_back(std::make_pair(dataratio, "histl"));
    
    cout << "Mults" << endl;


    // ----------------------------------------------------------------------------------
    //
    // First make pt in bins of pz

    // Values to multiply each bin by to get them on a similar range
    //These are for the XY!! 
    vector<double> multipliers1 = GetScales(histAndOpts,false,4.59,0.75,false);
    bool doMultipliers = false;
    GridCanvas* gc=NULL;
    gc=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV)", "#Sigma T_{p}", doMultipliers ? &multipliers1[0] : NULL);
    // Set the y range manually. Can also use gc->Remax() to guess automatically
    gc->SetYLimits(0.01,2.49);
    gc->SetYTitle(Form("Ratio to %s",sample2.c_str()));
    gc->Modified();
    // Example of adding a legend. The co-ordinate system is NDC on the
    // entire canvas, ie (0,0) in the bottom left corner of the canvas
    // (not the individual pad), and (1,1) in the top right
    TLegend* leg=new TLegend(0.8, 0.1, 1, 0.3);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(data1results[i], sample1.c_str(), "lpe");
    leg->AddEntry(data2results[i], sample2.c_str(), "lpe");
    

    TLatex mytex;
    mytex.SetTextSize(0.05);
    string mystring =     Form("%.2f < P_{||} [GeV] < %.2f",pz3Dbins[i-1],pz3Dbins[i]);
    mytex.DrawLatex(0.35,0.96,mystring.c_str());
    gc->Print(doMultipliers ? Form("nu-3d-datacomp-xsec-comps-pt-multiplier_bin_%d.eps",i) : Form("nu-3d-datacomp-xsec-comps-pt_bin_%d.eps",i));
    gc->Print(doMultipliers ? Form("nu-3d-datacomp-xsec-comps-pt-multiplier_bin_%d.png",i) : Form("nu-3d-datacomp-xsec-comps-pt_bin_%d.png",i));
    gc->Print(doMultipliers ? Form("nu-3d-datacomp-xsec-comps-pt-multiplier_bin_%d.C",i) : Form("nu-3d-datacomp-xsec-comps-pt_bin_%d.C",i));
      
    // // ----------------------------------------------------------------------------------
    // //
    // // Now make pz in bins of pt. It's all the same

    // // Values to multiply each bin by to get them on a similar range
    vector<double> multiplierspz1 = GetScales(histAndOpts,true,5.49,0.75,false);
    
    // // plotXAxis1D fiddles the x axis values to squash up the tail so it
    // // doesn't take up all the horizontal space.
    GridCanvas* gc2= NULL;
    gc2=plotXAxis1D(histAndOpts, "#Sigma T_{p} (GeV)", "P_{t,muon}", doMultipliers ? &multiplierspz1[0] : NULL);

    gc2->SetYLimits(0.01,2.49);
    gc2->SetYTitle(Form("Ratio to %s",sample2.c_str()));
    gc2->Modified();
    mytex.DrawLatex(0.35,0.96,mystring.c_str());
    gc2->Print(doMultipliers ? Form("nu-3d-datacomp-xsec-comps-pz-multiplier_bin_%d.eps",i) : Form("nu-3d-datacomp-xsec-comps-pz_bin_%d.eps",i));
    gc2->Print(doMultipliers ? Form("nu-3d-datacomp-xsec-comps-pz-multiplier_bin_%d.png",i) : Form("nu-3d-datacomp-xsec-comps-pz_bin_%d.png",i));
    gc2->Print(doMultipliers ? Form("nu-3d-datacomp-xsec-comps-pz-multiplier_bin_%d.C",i) : Form("nu-3d-datacomp-xsec-comps-pz_bin_%d.C",i));

    GridCanvas* gc3= NULL;
    gc3=plotXAxis1D(histAndOpts2, "#Sigma T_{p} (GeV)", "P_{t,muon}", doMultipliers ? &multiplierspz1[0] : NULL);

    gc3->SetYLimits(0.6,1.4);
    gc3->SetYTitle(Form("Ratio to %s",sample2.c_str()));
    gc3->Modified();
    mytex.DrawLatex(0.35,0.96,mystring.c_str());
    gc3->Print(doMultipliers ? Form("nu-3d-datacomp-xsec-comps-pz-multiplier_bin_%d_dataRatio.eps",i) : Form("nu-3d-datacomp-xsec-comps-pz_bin_%d_dataRatio.eps",i));
    gc3->Print(doMultipliers ? Form("nu-3d-datacomp-xsec-comps-pz-multiplier_bin_%d_dataRatio.png",i) : Form("nu-3d-datacomp-xsec-comps-pz_bin_%d_dataRatio.png",i));
    gc3->Print(doMultipliers ? Form("nu-3d-datacomp-xsec-comps-pz-multiplier_bin_%d_dataRatio.C",i) : Form("nu-3d-datacomp-xsec-comps-pz_bin_%d_dataRatio.C",i));
  
  }
}

int main(int argc, char* argv[])
{


  makePlots(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);

  return 0;
}
