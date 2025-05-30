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


void makePlots(bool doMultipliers,string location,string varx, string vary, bool doRatio, string sample="Standard", int projx_column=-1, int projx_row=-1, int projx_pixelx=-1, int projx_pixely=-1, int projy_column=-1, int projy_row=-1, int projy_pixelx=-1, int projy_pixely=-1)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  string var = varx+"_"+vary;
  string varspec = var;
  cout << "Variable Combination " << varspec << endl;
  string filename = Form("%s",location.c_str());//Final result
  TFile f1(filename.c_str());

  MnvH2D* dataMnv=(MnvH2D*)f1.Get(Form("h_%s_data_nobck_unfold_effcor_cross_section",varspec.c_str()));
  MnvH2D* mcMnv=(MnvH2D*)f1.Get(Form("h_%s_mc_nobck_unfold_effcor_cross_section",varspec.c_str()));
  MnvH2D* mcMnv_qelike_qe = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_qe",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_res = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_res",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_dis = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_dis",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_2p2h = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_2p2h",varspec.c_str()));//Get from N track

  MnvH2D* mcMnv_qelike_qe_protonfsi = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_qe_protonfsi",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_qe_neutronfsi = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_qe_neutronfsi",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_res_protonfsi = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_res_protonfsi",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_res_neutronfsi = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_res_neutronfsi",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_res_protonisi = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_res_protonisi",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_res_neutronisi = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_res_neutronisi",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_2p2h_nn = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_2p2h_nn",varspec.c_str()));//Get from N track
  MnvH2D* mcMnv_qelike_2p2h_np = (MnvH2D*)f1.Get(Form("h_%s_cross_section_qelike_2p2h_np",varspec.c_str()));//Get from N track

  cout << "Loaded " << endl;

  dataMnv->Scale(1e39, "width");
  mcMnv->Scale(1e39, "width");
  cout << "A" << endl;
  mcMnv_qelike_qe->Scale(1e39,"width");
  cout << "A1" << endl;
  mcMnv_qelike_qe_protonfsi->Scale(1e39,"width");
  cout << "A2" << endl;
  mcMnv_qelike_qe_neutronfsi->Scale(1e39,"width");
  cout << "B" << endl;
  mcMnv_qelike_res->Scale(1e39,"width");
  mcMnv_qelike_res_protonfsi->Scale(1e39,"width");
  mcMnv_qelike_res_neutronfsi->Scale(1e39,"width");
  mcMnv_qelike_res_protonisi->Scale(1e39,"width");
  mcMnv_qelike_res_neutronisi->Scale(1e39,"width");
  mcMnv_qelike_dis->Scale(1e39,"width");
  mcMnv_qelike_2p2h->Scale(1e39,"width");
  mcMnv_qelike_2p2h_nn->Scale(1e39,"width");
  mcMnv_qelike_2p2h_np->Scale(1e39,"width");
  
  cout << "Scaled" << endl;

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());

  TH2* mc_qelike_qe = new TH2D(mcMnv_qelike_qe->GetCVHistoWithStatError());
  TH2* mc_qelike_qe_protonfsi = new TH2D(mcMnv_qelike_qe_protonfsi->GetCVHistoWithStatError());
  TH2* mc_qelike_qe_neutronfsi = new TH2D(mcMnv_qelike_qe_neutronfsi->GetCVHistoWithStatError());
  TH2* mc_qelike_res = new TH2D(mcMnv_qelike_res->GetCVHistoWithStatError());
  TH2* mc_qelike_res_protonfsi = new TH2D(mcMnv_qelike_res_protonfsi->GetCVHistoWithStatError());
  TH2* mc_qelike_res_neutronfsi = new TH2D(mcMnv_qelike_res_neutronfsi->GetCVHistoWithStatError());
  TH2* mc_qelike_res_protonisi = new TH2D(mcMnv_qelike_res_protonisi->GetCVHistoWithStatError());
  TH2* mc_qelike_res_neutronisi = new TH2D(mcMnv_qelike_res_neutronisi->GetCVHistoWithStatError());
  TH2* mc_qelike_dis = new TH2D(mcMnv_qelike_dis->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h = new TH2D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h_nn = new TH2D(mcMnv_qelike_2p2h_nn->GetCVHistoWithStatError());
  TH2* mc_qelike_2p2h_np = new TH2D(mcMnv_qelike_2p2h_np->GetCVHistoWithStatError());

  cout << "Ratio? " <<endl;
  

  if(doRatio){
    dataStat->Divide(mc);
    data->Divide(mc);
    mc_qelike_qe->Divide(mc);
    mc_qelike_qe_protonfsi->Divide(mc);
    mc_qelike_qe_neutronfsi->Divide(mc);
    mc_qelike_res->Divide(mc);
    mc_qelike_res_protonfsi->Divide(mc);
    mc_qelike_res_neutronfsi->Divide(mc);
    mc_qelike_res_protonisi->Divide(mc);
    mc_qelike_res_neutronisi->Divide(mc);
    mc_qelike_dis->Divide(mc);
    mc_qelike_2p2h->Divide(mc);
    mc_qelike_2p2h_nn->Divide(mc);
    mc_qelike_2p2h_np->Divide(mc);
    mc->Divide(mc);
  }


  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);
  

  mc_qelike_qe->SetLineColor(mycolors[3]);
  mc_qelike_qe_protonfsi->SetLineColor(mycolors[3]);
  mc_qelike_qe_neutronfsi->SetLineColor(mycolors[3]);
  mc_qelike_res->SetLineColor(mycolors[4]);
  mc_qelike_res_protonfsi->SetLineColor(mycolors[4]);
  mc_qelike_res_neutronfsi->SetLineColor(mycolors[4]);
  mc_qelike_res_protonisi->SetLineColor(mycolors[4]);
  mc_qelike_res_neutronisi->SetLineColor(mycolors[4]);
  mc_qelike_dis->SetLineColor(mycolors[5]);
  mc_qelike_2p2h->SetLineColor(mycolors[16]);
  mc_qelike_2p2h_nn->SetLineColor(mycolors[16]);
  mc_qelike_2p2h_np->SetLineColor(mycolors[16]);


  mc_qelike_qe_protonfsi->SetLineStyle(2);
  mc_qelike_qe_neutronfsi->SetLineStyle(3);

  mc_qelike_res_protonfsi->SetLineStyle(2);
  mc_qelike_res_neutronfsi->SetLineStyle(3);
  mc_qelike_res_protonisi->SetLineStyle(2);
  mc_qelike_res_neutronisi->SetLineStyle(3);

  mc_qelike_2p2h_nn->SetLineStyle(2);
  mc_qelike_2p2h_np->SetLineStyle(3);


  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(1);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);
  dataStat->SetLineColor(kBlack);
  dataStat->SetMarkerStyle(1);

  // Make a list of the histograms we want to draw, along with the
  // draw options we want to use for them. You can add "graph" to the
  // draw options if you want the histogram to be converted to a graph
  // and then drawn. In that case the draw options are interpreted as
  // options to TGraphErrors::Draw().
  //
  // I don't know what happens if you put a "graph" first in the list,
  // so don't do that. Make sure the first item doesn't have "graph"
  // in its options
  cout << "Setup histogram vectors " << sample << endl;
  std::vector<std::pair<TH2*, const char*> > histAndOpts;
  histAndOpts.push_back(std::make_pair(dataStat, "histpe1"));
  histAndOpts.push_back(std::make_pair(mc,       "hist l"));
  if(sample=="Standard"){
    histAndOpts.push_back(std::make_pair(mc_qelike_qe,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_dis,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h,       "hist l"));
  }
  else if(sample=="QE"){
    histAndOpts.push_back(std::make_pair(mc_qelike_qe,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_qe_protonfsi,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_qe_neutronfsi,       "hist l"));
  }
  else if(sample=="RESFSI"){
    histAndOpts.push_back(std::make_pair(mc_qelike_res,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res_protonfsi,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res_neutronfsi,       "hist l"));
  }
  else if(sample=="RESISI"){
    histAndOpts.push_back(std::make_pair(mc_qelike_res,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res_protonisi,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res_neutronisi,       "hist l"));
  }
 
  else if(sample=="2p2h"){
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_nn,       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_np,       "hist l"));
  }
  else exit(1);

  histAndOpts.push_back(std::make_pair(data,     "histpe1"));


  cout << "Done Setting up vectors" << endl;

  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  vector<double> multi_1 = GetScales(histAndOpts, false,4.99,0.75);

  GridCanvas* gc=NULL;
  if(projx_column==-1){
    gc= plotYAxis1D(histAndOpts, vary.c_str() ,varx.c_str() ,doMultipliers ? &multi_1[0] : NULL);
  }
  else{
    gc= plotYAxis1D(histAndOpts, vary.c_str() ,varx.c_str(),projy_column,projy_row,projy_pixelx,projy_pixely ,doMultipliers ? &multi_1[0] : NULL);
  }
  
  if(doRatio){
    gc->SetYLimits(0, 1.99);
    gc->SetYTitle("Ratio to Minerva Tune v1.0.1");
  }
  else{
    gc->SetYLimits(0, 4.99);
    gc->SetYTitle("d^{2}#sigma/d (x10^{-39} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  }
    
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  //  TLegend* leg=new TLegend(0.17, 0.7, 0.31, 0.9);
  TLegend* leg=new TLegend(0.7, 0.1, 0.9, 0.3);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MnvGENIE", "l");
  if(sample=="Standard"){
    leg->AddEntry(mc_qelike_qe,"QE","l");
    leg->AddEntry(mc_qelike_res,"Resonant","l");
    leg->AddEntry(mc_qelike_dis,"DIS","l");
    leg->AddEntry(mc_qelike_2p2h,"2p2h","l");
  }
  else if(sample=="QE"){
    leg->AddEntry(mc_qelike_qe,"QE","l");
    leg->AddEntry(mc_qelike_qe_protonfsi,"QE Proton Dominant","l");
    leg->AddEntry(mc_qelike_qe_neutronfsi,"QE Neutron Dominant","l");
  }    
  else if(sample=="RESFSI"){
    leg->AddEntry(mc_qelike_res,"Resonant","l");
    leg->AddEntry(mc_qelike_res_protonfsi,"Resonant Proton Dominant","l");
    leg->AddEntry(mc_qelike_res_neutronfsi,"Resonant Neutron Dominant","l");
  }
  else if(sample=="RESISI"){
    leg->AddEntry(mc_qelike_res,"Resonant","l");
    leg->AddEntry(mc_qelike_res_protonisi,"Resonant Proton Dominant","l");
    leg->AddEntry(mc_qelike_res_neutronisi,"Resonant Neutron Dominant","l");
  }
  else if(sample=="2p2h"){
    leg->AddEntry(mc_qelike_2p2h,"2p2h","l");
    leg->AddEntry(mc_qelike_2p2h_nn,"2p2h nn-pair","l");
    leg->AddEntry(mc_qelike_2p2h_np,"2p2h np-pair","l");
  }

  

  if(doRatio){
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-multiplier_ratio.eps",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s_ratio.eps",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-multiplier_ratio.png",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s_ratio.png",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-multiplier_ratio.C",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s_ratio.C",var.c_str()));
  }
  else{
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s.eps",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s.png",var.c_str()));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-projy-%s-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-projy-%s.C",var.c_str()));
  }
  

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  vector<double> multi_2 = GetScales(histAndOpts,true,4.99,0.75);
  


  // Values to multiply each bin by to get them on a similar range


  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  GridCanvas* gc2=NULL;
  if(projx_column==-1){
    gc2 =plotXAxis1D(histAndOpts, varx.c_str() , vary.c_str() ,doMultipliers ? &multi_2[0] : NULL);
  }
  else{
    gc2 =plotXAxis1D(histAndOpts, varx.c_str() , vary.c_str(), projx_column, projx_row, projx_pixelx, projx_pixely, doMultipliers ? &multi_2[0] : NULL);
  }
  if(doRatio) {
    gc2->SetYLimits(0, 1.99);
    gc2->SetYTitle("Ratio to Minerva Tune v1.0.1");

  }
  else{
    gc2->SetYLimits(0, 4.99);
    gc2->SetYTitle("d^{2}#sigma/d_{t} (x10^{-39} cm^{2}/GeV^{2}/c^{2}/Nucleon)");
  }


  gc2->Modified();
  if(doRatio){
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-multiplier_ratio.eps",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s_ratio.eps",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-multiplier_ratio.png",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s_ratio.png",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-multiplier_ratio.C",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s_ratio.C",var.c_str()));
  }
  else{
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-multiplier.eps",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s.eps",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s.png",var.c_str()));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-projx-%s-multiplier.C",var.c_str()) : Form("nu-2d-xsec-comps-projx-%s.C",var.c_str()));

  }
  

}

int main(int argc, char* argv[])
{
  if(argc==5){
    makePlots(true,argv[1],argv[2],argv[3],false,argv[4]);
    makePlots(false,argv[1],argv[2],argv[3],true,argv[4]);
  }
  if(argc>5){
    makePlots(true,argv[1],argv[2],argv[3],false,argv[4],atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9]),atoi(argv[10]),atoi(argv[11]),atoi(argv[12]));
    makePlots(false,argv[1],argv[2],argv[3],true,argv[4],atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9]),atoi(argv[10]),atoi(argv[11]),atoi(argv[12]));
  }

  return 0;
}
