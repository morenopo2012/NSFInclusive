//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/HyperDimLinearizer.h"//THIS HAS TO CHANGE TO BE INCLUDED IN THE MAKE FILE EVENTUALLY.
#include "PlotUtils/MnvColors.h"
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
#include "TLine.h"
#include "TArrow.h"
//#include "localColor.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;



void makePlots(bool doMultipliers,bool doRatio, string location, string location2, string location3, string location4, string location5)
{
  
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  string variableset = "pzptrec";
  string yaxislabel = "P_{t} (GeV/c)";
  string xaxislabel = "P_{||} (GeV/c)";
  string yaxislong = "Muon Transverse Momentum (GeV/c)";
  string crosssectionlabel = "d^{3}#sigma/dp_{T}dp_{||}d#SigmaT_{p} (x10^{-39} cm^{2}/GeV^{3}/c^{2}/Nucleon)";
  
  string histname_data = "xsec_data";
  string histname_mc = "xsec_mc";
  
  TFile f_ch(Form("%s",location.c_str()));//final data
  TFile f_pb(Form("%s",location2.c_str()));//final mc
  TFile f_fe(Form("%s",location3.c_str()));//final mc
  TFile f_c(Form("%s",location4.c_str()));//final mc
  TFile f_h20(Form("%s",location5.c_str()));//final mc
  
  //ch
  MnvH2D* dataMnv_ch=(MnvH2D*)f_ch.Get(histname_data.c_str());
  MnvH2D* mcMnv_ch=(MnvH2D*)f_ch.Get(histname_mc.c_str());
  //pb
  MnvH2D* dataMnv_pb=(MnvH2D*)f_pb.Get(histname_data.c_str());
  MnvH2D* mcMnv_pb=(MnvH2D*)f_pb.Get(histname_mc.c_str());
  //fe
  MnvH2D* dataMnv_fe=(MnvH2D*)f_fe.Get(histname_data.c_str());
  MnvH2D* mcMnv_fe=(MnvH2D*)f_fe.Get(histname_mc.c_str());
  //c
  MnvH2D* dataMnv_c=(MnvH2D*)f_c.Get(histname_data.c_str());
  MnvH2D* mcMnv_c=(MnvH2D*)f_c.Get(histname_mc.c_str());
  //h20
  MnvH2D* dataMnv_h20=(MnvH2D*)f_h20.Get(histname_data.c_str());
  MnvH2D* mcMnv_h20=(MnvH2D*)f_h20.Get(histname_mc.c_str());
  
  
  int n_material = 5;
  
  std::vector<TH2D*> data_materials;
  std::vector<TH2D*> mc_materials;






  data_materials.push_back(new TH2D(dataMnv_h20->GetCVHistoWithError()));
  data_materials.push_back(new TH2D(dataMnv_c->GetCVHistoWithError()));
  data_materials.push_back(new TH2D(dataMnv_fe->GetCVHistoWithError()));
  data_materials.push_back(new TH2D(dataMnv_pb->GetCVHistoWithError()));
  data_materials.push_back(new TH2D(dataMnv_ch->GetCVHistoWithError()));


  mc_materials.push_back(new TH2D(mcMnv_h20->GetCVHistoWithStatError()));
  mc_materials.push_back(new TH2D(mcMnv_c->GetCVHistoWithStatError()));
  mc_materials.push_back(new TH2D(mcMnv_fe->GetCVHistoWithStatError()));
  mc_materials.push_back(new TH2D(mcMnv_pb->GetCVHistoWithStatError()));
  mc_materials.push_back(new TH2D(mcMnv_ch->GetCVHistoWithStatError()));

  vector<string> material;




  material.push_back("H_{2}O");
  material.push_back("C");
  material.push_back("Fe");
  material.push_back("Pb");
  material.push_back("CH");


  bool dooverall = true;
  vector<double> overallscale;
  if(!doMultipliers){
    overallscale.push_back(1);
    overallscale.push_back(1);
    overallscale.push_back(1);
    overallscale.push_back(1);
    overallscale.push_back(1);
    overallscale.push_back(1);
  }
  else{
    for(unsigned int i=0;i<n_material;i++) overallscale.push_back(1);
  }
  
  for(unsigned int i=0;i<n_material;i++){
    data_materials[i]->Scale(1e39*overallscale[i],"width");
    mc_materials[i]->Scale(1e39*overallscale[i],"width");
  }

    


  vector<int> mycolors = MnvColors::GetColors(10);//Light
  vector<int> mycolors2 = MnvColors::GetColors(11);//Dark
  vector<int> index;
  index.push_back(0);
  index.push_back(1);
  index.push_back(2);
  index.push_back(4);
  index.push_back(5);
  index.push_back(6);
  index.push_back(7);
  index.push_back(8);
  index.push_back(9);
  index.push_back(10);
  
  std::vector<std::pair<TH2*, const char*> > histAndOpts;
  std::vector<std::pair<TH2*, const char*> > histAndOptsDataPoints;
  TLegend *leg = new TLegend(0.15,0.95,0.9,1);//if I do 3 separate legends
  leg->SetNColumns(n_material);
  vector<double> scales ;  
  for(unsigned int i=0;i<n_material-1;i++){//skip underflow and overflow pz bins
    cout << "DOING BIN " << i << endl;
    // Get the data histogram with stat error and with total error
    // separately so we can plot them both for inner and outer ticks
    // These line and marker styles will be propagated to the 1D plots
    mc_materials[i]->SetLineColor(mycolors2[index[i]]);
    mc_materials[i]->SetLineWidth(2);
    leg->AddEntry(data_materials[i],material[i].c_str(),"P");
    
    data_materials[i]->SetLineColor(mycolors[index[i]]);
    //    data_materials[i]->SetLineWidth(i+1);
    data_materials[i]->SetMarkerColor(mycolors[index[i]]);
    data_materials[i]->SetMarkerSize(0.7);
    
    data_materials[i]->SetFillColor(mycolors[index[i]]);
    data_materials[i]->SetFillStyle(3003);
    
    
    if(doRatio){
      TF1 *myconstant = new TF1("myfunc","1",0,50);
      data_materials[i]->Divide(data_materials[4]);
      mc_materials[i]->Divide(mc_materials[4]);

      
      double sepsize=1;
      data_materials[i]->Add(myconstant,(i)*sepsize);
      mc_materials[i]->Add(myconstant,(i)*sepsize);
      
      
    }

    //      data_materials_statonly_simple[i]->SetMarkerColor(mycolors2[index[i]]);
    // Make a list of the histograms we want to draw, along with the
    // draw options we want to use for them. You can add "graph" to the
    // draw options if you want the histogram to be converted to a graph
    // and then drawn. In that case the draw options are interpreted as
    // options to TGraphErrors::Draw().
    //
    // I don't know what happens if you put a "graph" first in the list,
    // so don't do that. Make sure the first item doesn't have "graph"
    // in its options
    
    
    histAndOpts.push_back(std::make_pair(data_materials[i], "histpe1"));
    histAndOpts.push_back(std::make_pair(mc_materials[i],       "histl"));
    /*
    if(!doRatio){
      histAndOpts.push_back(std::make_pair(data_materials[i], "histle5"));
    }
    else{
      histAndOpts.push_back(std::make_pair(data_materials[i], "histle5"));
    }
    */
  }
  cout << "Mults" << endl;
  
  histAndOpts.insert(histAndOpts.end(),histAndOptsDataPoints.begin(),histAndOptsDataPoints.end());
  
  // // ----------------------------------------------------------------------------------
  // //
  // // Now make pz in bins of pt. It's all the same
  
  // // Values to multiply each bin by to get them on a similar range
  
  // // plotXAxis1D fiddles the x axis values to squash up the tail so it
  // // doesn't take up all the horizontal space.
  scales = GetScales(histAndOpts,false,5.49,0.75,false);
  GridCanvas* gc= NULL;
  gc=plotYAxis1D_ReducedXRange(histAndOpts, yaxislabel, xaxislabel,4,3,0,2.2,800,500, doRatio?NULL:&scales[0]);
  //gc=plotYAxis1D(histAndOpts, yaxislabel, xaxislabel,doRatio?NULL:&scales[0]);
  gc->SetYTitle(crosssectionlabel.c_str());
  if(!doMultipliers){
    if(!doRatio) gc->SetYLimits(0,5.49);
    else{
      gc->SetYTitle("Ratio to A to CH result");
      gc->SetYLimits(0,6.99);
      //      gc->SetManualYLabels(nLabels, positions, valueStrings,0.03);
    }
  }
  else{
    gc->SetYLimits(0.00003, 10.49);
  }
  
  //  if(!doRatio)gc->SetLogy(true);
  leg->Draw("SAME");
  gc->Modified();
  if(!doRatio){
    gc->Print(doMultipliers ? Form("nu-nuclear-target-materials-%s-pz-multiplier_bin_log_all3D.eps",variableset.c_str()) : Form("nu-nuclear-target-materials-%s-pz_bin_log_all3D.eps",variableset.c_str()));
    gc->Print(doMultipliers ? Form("nu-nuclear-target-materials-%s-pz-multiplier_bin_log_all3D.png",variableset.c_str()) : Form("nu-nuclear-target-materials-%s-pz_bin_log_all3D.png",variableset.c_str()));
    gc->Print(doMultipliers ? Form("nu-nuclear-target-materials-%s-pz-multiplier_bin_log_all3D.C",variableset.c_str()) : Form("nu-nuclear-target-materials-%s-pz_bin_log_all3D.C",variableset.c_str()));
  }
  else{
    gc->Print(doMultipliers ? Form("nu-nuclear-target-materials-%s-pz-multiplier_bin_log_all3D_ratio.eps",variableset.c_str()) : Form("nu-nuclear-target-materials-%s-pz_bin_log_all3D_ratio.eps",variableset.c_str()));
    gc->Print(doMultipliers ? Form("nu-nuclear-target-materials-%s-pz-multiplier_bin_log_all3D_ratio.png",variableset.c_str()) : Form("nu-nuclear-target-materials-%s-pz_bin_log_all3D_ratio.png",variableset.c_str()));
    gc->Print(doMultipliers ? Form("nu-nuclear-target-materials-%s-pz-multiplier_bin_log_all3D_ratio.C",variableset.c_str()) : Form("nu-nuclear-target-materials-%s-pz_bin_log_all3D_ratio.C",variableset.c_str()));
    
  }
}

int main(int argc, char* argv[])
{
  
  //no ratio
  makePlots(true,false,argv[1],argv[2],argv[3],argv[4],argv[5]);
  makePlots(false,false,argv[1],argv[2],argv[3],argv[4],argv[5]);
  
  //ratio
  makePlots(false,true,argv[1],argv[2],argv[3],argv[4],argv[5]);

  return 0;
}
