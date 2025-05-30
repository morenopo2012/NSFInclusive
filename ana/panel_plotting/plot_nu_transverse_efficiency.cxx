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
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

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


void makePlots(bool doMultipliers,int anatype,string location,string varx, string vary)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  string var = varx+"_"+vary;
  string varspecfilename = var;
  string varspec = var;
  string filename = Form("%s_CV_Elastic_FSI/efficiency_combined/EffPurity_MakeFlux-1_%s_CombinedPlaylists.root",location.c_str(),varspec.c_str());
  TFile f1(filename.c_str());
  MnvH2D* mcMnv=(MnvH2D*)f1.Get(Form("h_%s_truth_qelike",varspec.c_str()));
  MnvH2D* mcMnv_qelike_qe = (MnvH2D*)f1.Get(Form("h_%s_truth_qelike_qe",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_res = (MnvH2D*)f1.Get(Form("h_%s_truth_qelike_res",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_dis = (MnvH2D*)f1.Get(Form("h_%s_truth_qelike_dis",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_2p2h = (MnvH2D*)f1.Get(Form("h_%s_truth_qelike_2p2h",varspec.c_str()));//Get from N track

  MnvH2D* mcMnv_qelike_qe_proton_fsi = (MnvH2D*)f1.Get(Form("h_%s_truth_qelike_qe_proton_fsi",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_qe_neutron_fsi = (MnvH2D*)f1.Get(Form("h_%s_truth_qelike_qe_neutron_fsi",varspec.c_str()));//Get from N track

  MnvH2D* mcMnv_qelike_res_proton_fsi = (MnvH2D*)f1.Get(Form("h_%s_truth_qelike_res_proton_fsi",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_res_neutron_fsi = (MnvH2D*)f1.Get(Form("h_%s_truth_qelike_res_neutron_fsi",varspec.c_str()));//Get from N track

  MnvH2D* mcMnv_qelike_2p2h_np = (MnvH2D*)f1.Get(Form("h_%s_truth_qelike_2p2h_np",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_2p2h_nn = (MnvH2D*)f1.Get(Form("h_%s_truth_qelike_2p2h_nn",varspec.c_str()));//Get from N track



  mcMnv->Scale(1e-5, "width");
  mcMnv_qelike_qe->Scale(1e-5,"width");
  mcMnv_qelike_res->Scale(1e-5,"width");
  mcMnv_qelike_dis->Scale(1e-5,"width");
  mcMnv_qelike_2p2h->Scale(1e-5,"width");
  mcMnv_qelike_qe_proton_fsi->Scale(1e-5,"width");
  mcMnv_qelike_qe_neutron_fsi->Scale(1e-5,"width");
  mcMnv_qelike_res_proton_fsi->Scale(1e-5,"width");
  mcMnv_qelike_res_neutron_fsi->Scale(1e-5,"width");
  mcMnv_qelike_2p2h_np->Scale(1e-5,"width");
  mcMnv_qelike_2p2h_nn->Scale(1e-5,"width");


  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mc_qelike_qe = new TH2D(mcMnv_qelike_qe->GetCVHistoWithStatError());
  TH2* mc_qelike_res = new TH2D(mcMnv_qelike_res->GetCVHistoWithStatError());
  TH2* mc_qelike_dis = new TH2D(mcMnv_qelike_dis->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h = new TH2D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());

  TH2* mc_qelike_qe_proton_fsi = new TH2D(mcMnv_qelike_qe_proton_fsi->GetCVHistoWithStatError());
  TH2* mc_qelike_qe_neutron_fsi = new TH2D(mcMnv_qelike_qe_neutron_fsi->GetCVHistoWithStatError());

  TH2* mc_qelike_res_proton_fsi = new TH2D(mcMnv_qelike_res_proton_fsi->GetCVHistoWithStatError());
  TH2* mc_qelike_res_neutron_fsi = new TH2D(mcMnv_qelike_res_neutron_fsi->GetCVHistoWithStatError());

  TH2* mc_qelike_2p2h_np = new TH2D(mcMnv_qelike_2p2h_np->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h_nn = new TH2D(mcMnv_qelike_2p2h_nn->GetCVHistoWithStatError());



  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);

  mc_qelike_qe->SetLineColor(mycolors[3]);
  mc_qelike_res->SetLineColor(mycolors[4]);
  mc_qelike_dis->SetLineColor(mycolors[5]);
  mc_qelike_2p2h->SetLineColor(mycolors[16]);

  mc_qelike_qe_proton_fsi->SetLineColor(mycolors[3]);
  mc_qelike_qe_neutron_fsi->SetLineColor(mycolors[3]);
  mc_qelike_qe_proton_fsi->SetLineStyle(2);
  mc_qelike_qe_neutron_fsi->SetLineStyle(3);

  mc_qelike_res_proton_fsi->SetLineColor(mycolors[4]);
  mc_qelike_res_neutron_fsi->SetLineColor(mycolors[4]);
  mc_qelike_res_proton_fsi->SetLineStyle(2);
  mc_qelike_res_neutron_fsi->SetLineStyle(3);

  mc_qelike_2p2h_np->SetLineColor(mycolors[16]);
  mc_qelike_2p2h_nn->SetLineColor(mycolors[16]);
  mc_qelike_2p2h_np->SetLineStyle(2);
  mc_qelike_2p2h_nn->SetLineStyle(3);


  // These line and marker styles will be propagated to the 1D plots


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
  histAndOpts.push_back(std::make_pair(mc,       "hist l"));
  if(anatype==1){
    histAndOpts.push_back(std::make_pair(mc_qelike_qe,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_dis,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h,       "hist l"));
  }
  else if(anatype==2){
    histAndOpts.push_back(std::make_pair(mc_qelike_qe,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_qe_proton_fsi,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_qe_neutron_fsi,       "hist l"));


  }
  else if(anatype==3){
    histAndOpts.push_back(std::make_pair(mc_qelike_res,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res_proton_fsi,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res_neutron_fsi,       "hist l"));
  }
  else if(anatype==4){
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_np,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_nn,       "hist l"));
  }
  
  



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
  gc->SetYTitle("Event Rate (scaled by 1e-3)");
    
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  //  TLegend* leg=new TLegend(0.17, 0.7, 0.31, 0.9);
  TLegend* leg=new TLegend(0.7, 0.1, 0.9, 0.3);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(mc, "MnvGENIE v1.0.1", "l");
  leg->AddEntry(mc_qelike_qe,"QE","l");
  leg->AddEntry(mc_qelike_res,"Resonant","l");
  leg->AddEntry(mc_qelike_dis,"DIS","l");
  leg->AddEntry(mc_qelike_2p2h,"2p2h","l");
  
  
  if(anatype==1){
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s.eps",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s.png",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s.C",var.c_str()));
  }
  else if(anatype==2){
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-qetype-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s-qetype.eps",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-qetype-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s-qetype.png",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-qetype-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s-qetype.C",var.c_str()));

  }
  
  else if(anatype==3){
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-restype-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s-restype.eps",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-restype-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s-restype.png",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-restype-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s-restype.C",var.c_str()));

  }
  else if(anatype==4){
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-2p2htype-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s-2p2htype.eps",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-2p2htype-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s-2p2htype.png",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-2p2htype-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s-2p2htype.C",var.c_str()));

  }

  

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
  gc2->SetYTitle("Event Rate (scaled by 1e-3)");

  gc2->Modified();

  if(anatype==1){
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s.eps",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s.png",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s.C",var.c_str()));
  }
  else if(anatype==2){
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-qetype-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s-qetype.eps",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-qetype-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s-qetype.png",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-qetype-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s-qetype.C",var.c_str()));

  }
  else if(anatype==3){
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-restype-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s-restype.eps",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-restype-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s-restype.png",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-restype-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s-restype.C",var.c_str()));
    
  }
  
  
  else if(anatype==4){
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-2p2htype-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s-2p2htype.eps",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-2p2htype-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s-2p2htype.png",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-2p2htype-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s-2p2htype.C",var.c_str()));
    
  }  
  }
  

int main(int argc, char* argv[])
{
  // makePlots(true,true,argv[1],argv[2]);
  // makePlots(true,false,argv[1],argv[2]);
  // makePlots(false,true,argv[1],argv[2]);
  //  makePlots(false,false,argv[1],argv[2],argv[3]);
  makePlots(true,1,argv[1],argv[2],argv[3]);//standard
  makePlots(true,2,argv[1],argv[2],argv[3]);//qe HE particle
  makePlots(true,3,argv[1],argv[2],argv[3]);//res HE particle
  makePlots(true,4,argv[1],argv[2],argv[3]);//2p2h nn or np

  return 0;
}
