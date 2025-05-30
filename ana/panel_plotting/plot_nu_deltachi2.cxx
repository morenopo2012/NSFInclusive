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

void makePlots(string location)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  gStyle->SetLabelSize(0.04);

  TFile f1(Form("%s_CV/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");
  MnvH2D* mcMnv=(MnvH2D*)f1.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");

  TFile f2(Form("%s_deltaChi2/DeltaChi2ByBin.root",location.c_str()));


  dataMnv->GetXaxis()->SetLabelSize(0.04);
  dataMnv->GetYaxis()->SetLabelSize(0.04);


  MnvH2D* mcMnv_qe = (MnvH2D*)f1.Get("h_pzmu_ptmu_cross_section_qe");//Get from N track
  MnvH2D* mcMnv_res = (MnvH2D*)f1.Get("h_pzmu_ptmu_cross_section_res");//Get from N track
  MnvH2D* mcMnv_dis = (MnvH2D*)f1.Get("h_pzmu_ptmu_cross_section_dis");//Get from N track
  MnvH2D* mcMnv_dis_dis = (MnvH2D*)f1.Get("h_pzmu_ptmu_cross_section_dis_dis");//Get from N track
  MnvH2D* mcMnv_dis_sis = (MnvH2D*)f1.Get("h_pzmu_ptmu_cross_section_dis_sis");//Get from N track
  MnvH2D* mcMnv_2p2h = (MnvH2D*)f1.Get("h_pzmu_ptmu_cross_section_2p2h");//Get from N track
  MnvH2D* mcMnv_oth = (MnvH2D*)f1.Get("h_pzmu_ptmu_cross_section_oth");//Get from N track

  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV/c)");

  dataMnv->Scale(1e39, "width");
  mcMnv->Scale(1e39, "width");

  mcMnv_qe->Scale(1e39,"width");
  mcMnv_res->Scale(1e39,"width");
  mcMnv_dis->Scale(1e39,"width");
  mcMnv_dis_dis->Scale(1e39,"width");
  mcMnv_dis_sis->Scale(1e39,"width");
  mcMnv_2p2h->Scale(1e39,"width");
  mcMnv_oth->Scale(1e39,"width");

  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());

  TH2* mc_qe = new TH2D(mcMnv_qe->GetCVHistoWithStatError());
  TH2* mc_res = new TH2D(mcMnv_res->GetCVHistoWithStatError());
  TH2* mc_dis = new TH2D(mcMnv_dis->GetCVHistoWithStatError());
  TH2* mc_dis_dis = new TH2D(mcMnv_dis_dis->GetCVHistoWithStatError());
  TH2* mc_dis_sis = new TH2D(mcMnv_dis_sis->GetCVHistoWithStatError());
  TH2* mc_2p2h = new TH2D(mcMnv_2p2h->GetCVHistoWithStatError());
  TH2* mc_oth = new TH2D(mcMnv_oth->GetCVHistoWithStatError());

  TH2* dc2_genie2126 = (TH2*)f2.Get("genie2126withNieves");
  TH2* dc2_disamu = (TH2*)f2.Get("DIS_AMU");
  TH2* dc2_disncteq = (TH2*)f2.Get("DIS_NCTEQ");
  TH2* dc2_disncteqnu = (TH2*)f2.Get("DIS_NCTEQNu");
  TH2* dc2_mk = (TH2*)f2.Get("MK");
  TH2* dc2_rparesminos = (TH2*)f2.Get("RPA_Res_MINOS");
  TH2* dc2_rparesnieves = (TH2*)f2.Get("RPA_Res_Nieves");
  TH2* dc2_mnvgeniev2 = (TH2*)f2.Get("MnvGENIEv2");
  TH2* dc2_pion2p2h = (TH2*)f2.Get("pion_2p2h");
  TH2* dc2_pionrpa = (TH2*)f2.Get("pion_rpa");
  TH2* dc2_piontune = (TH2*)f2.Get("piontune");
  TH2* dc2_gibuu = (TH2*)f2.Get("GiBUU");
  TH2* dc2_nuwrosf = (TH2*)f2.Get("NuWro_SF");
  TH2* dc2_nuwrolfg = (TH2*)f2.Get("NuWro_LFG");
  TH2* dc2_neutlfg = (TH2*)f2.Get("NEUT_LFG_105");
  TH2* dc2_neutsf = (TH2*)f2.Get("NEUT_SF_103");

  vector<TH2*> dc2_vec;
  dc2_vec.push_back(dc2_genie2126);
  dc2_vec.push_back(dc2_disamu);
  dc2_vec.push_back(dc2_disncteq);
  dc2_vec.push_back(dc2_disncteqnu);
  dc2_vec.push_back(dc2_mk);
  dc2_vec.push_back(dc2_rparesminos);
  dc2_vec.push_back(dc2_rparesnieves);
  dc2_vec.push_back(dc2_mnvgeniev2);
  dc2_vec.push_back(dc2_pion2p2h);
  dc2_vec.push_back(dc2_pionrpa);
  dc2_vec.push_back(dc2_piontune);
  dc2_vec.push_back(dc2_gibuu);
  dc2_vec.push_back(dc2_nuwrosf);
  dc2_vec.push_back(dc2_nuwrolfg);
  dc2_vec.push_back(dc2_neutlfg);
  dc2_vec.push_back(dc2_neutsf);
  

  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);

  mc_qe->SetLineColor(mycolors[3]);
  mc_res->SetLineColor(mycolors[5]);
  mc_dis_dis->SetLineColor(kViolet-3);
  mc_dis_sis->SetLineColor(mycolors[4]);
  mc_oth->SetLineColor(mycolors[8]);


  //Add 2p2h with qe to reduce number of categories
  mc_qe->Add(mc_2p2h);

  mc_qe->SetLineWidth(1.5);
  mc_res->SetLineWidth(1.5);
  mc_dis_dis->SetLineWidth(1.5);
  mc_dis_sis->SetLineWidth(1.5);
  mc_oth->SetLineWidth(1.5);

  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.5);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);

  dataStat->SetLineColor(kBlack);


  //    TH2 *tmpden = (TH2*)mc->Clone("tmpden");
  //    tmpden->Sumw2(false);
  data->Divide(mc);
  dataStat->Divide(mc);
  mc_qe->Divide(mc);
  mc_res->Divide(mc);
  mc_dis_dis->Divide(mc);
  mc_dis_sis->Divide(mc);
  mc_oth->Divide(mc);
  mc->Divide(mc);
  
    

  for(int v=0;v<dc2_vec.size();v++){  

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
  histAndOpts.push_back(std::make_pair(dc2_vec[v], "hist"));


  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range

  GridCanvas* gc=plotXAxis1DRebinPz(histAndOpts, "Muon Longitudinal Momentum (GeV/c)", "p_{t}", 4,4,800,500,NULL);
  double max_val = dc2_vec[v]->GetMaximum();
  double min_val = dc2_vec[v]->GetMinimum();
  double ran = max_val;
  if(abs(max_val)<abs(min_val)) ran = abs(min_val);
  ran*=1.2;
  // Set the y range manually. Can also use gc->Remax() to guess automatically
  gc->SetYLimits(-1*ran,ran);
  gc->SetYTitle("#Delta#chi^{2}");
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  gc->Print(Form("nu-2d-deltachi2-%s-pt.eps",dc2_vec[v]->GetTitle()));
  gc->Print(Form("nu-2d-deltachi2-%s-pt.png",dc2_vec[v]->GetTitle()));
  gc->Print(Form("nu-2d-deltachi2-%s-pt.C",dc2_vec[v]->GetTitle()));
  
  

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  
  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  //GridCanvas* gc2=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV/c)", "p_{||}",4,4,800,500,NULL);
  GridCanvas* gc2=plotYAxis1DRebinPt(histAndOpts, "Muon Transverse Momentum (GeV/c)", "p_{||}",4,4,800,500,NULL);
  gc2->SetYLimits(-1*ran,ran);
  gc2->SetYTitle("#Delta#chi^{2}");
  gc2->Modified();
  gc2->Print(Form("nu-2d-deltachi2-%s-pz.eps",dc2_vec[v]->GetTitle()));
  gc2->Print(Form("nu-2d-deltachi2-%s-pz.png",dc2_vec[v]->GetTitle()));
  gc2->Print(Form("nu-2d-deltachi2-%s-pz.C",dc2_vec[v]->GetTitle()));
  
  }
}

int main(int argc, char* argv[])
{
  makePlots(argv[1]);
  return 0;
}
