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
//#include "localColor.h"
#include "PlotUtils/MnvColors.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;

void makePlots(bool doMultipliers,bool doRatio,string location, int modset=0, int areanorm=0)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  TFile f1(Form("%s_CV/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f2(Form("%s_default/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f3(Form("%s_CV_DIS_AMU/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f4(Form("%s_CV_DIS_NCTEQ/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f5(Form("%s_CV_MK/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f6(Form("%s_CV_RPA_Res_MINOS/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f7(Form("%s_CV_RPA_Res_Nieves/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f8(Form("%s_CV_minerva_joint_lowq2/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f9(Form("%s_pion_2p2h/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f10(Form("%s_pion_rpa/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f11(Form("%s_piontune/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f12(Form("%s_gibuu/GiBUU_MEInclusive2D_NotWidthScaled.root",location.c_str()));//Final result
  TFile f13(Form("%s_CV_DIS_NCTEQNU/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_pzmu.root",location.c_str()));//Final result
  TFile f14(Form("%s_NuWro/NuWro_SF.root",location.c_str()));//Final result
  TFile f15(Form("%s_NuWro/NuWro_LFG_bypart.root",location.c_str()));//Final result
  TFile f16(Form("%s_NEUT/NEUT_LFG_MA105.root",location.c_str()));//Final result
  TFile f17(Form("%s_NEUT/NEUT_SF_MA103.root",location.c_str()));//Final result
  TFile f18(Form("%s_gibuu_v2021/GiBUU_MEInclusive2D_NotWidthScaled_v2021.root",location.c_str()));//Final result

  //  TFile f18(Form("%s_NEUT/NEUT_SF_MA121.root",location.c_str()));//Final result
  //  TFile f14(Form("%s_NuWro/Tejin_NuWro.root",location.c_str()));//Final result

  MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section");
  MnvH2D* mcMnv=(MnvH2D*)f1.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");

  //Reweighted or different models
  MnvH2D* nommcMnv=(MnvH2D*)f2.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* amuMnv=(MnvH2D*)f3.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* ncteqMnv=(MnvH2D*)f4.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* ncteqnuMnv=(MnvH2D*)f13.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* mkMnv=(MnvH2D*)f5.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* resminosMnv=(MnvH2D*)f6.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* resnievesMnv=(MnvH2D*)f7.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* resminervaMnv=(MnvH2D*)f8.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* pion2p2hMnv=(MnvH2D*)f9.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* pionrpaMnv=(MnvH2D*)f10.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* piontuneMnv=(MnvH2D*)f11.Get("h_pzmu_ptmu_mc_nobck_unfold_effcor_cross_section");
  MnvH2D* gibuuMnv=new MnvH2D(*(TH2D*)f12.Get("GiBUU"));
  MnvH2D* gibuu_v2021Mnv=new MnvH2D(*(TH2D*)f18.Get("GiBUU"));
  MnvH2D* nuwrosfMnv= new MnvH2D(*(TH2D*)f14.Get("nuwro_SF"));
  MnvH2D* nuwrolfgMnv= new MnvH2D(*(TH2D*)f15.Get("nuwro_LFG"));
  MnvH2D* neutlfg105Mnv = new MnvH2D(*(TH2D*)f16.Get("ptpz"));
  MnvH2D* neutsf103Mnv = new MnvH2D(*(TH2D*)f17.Get("ptpz"));
  //  MnvH2D* neutsf121Mnv = new MnvH2D(*(TH2D*)f18.Get("ptpz"));



  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV/c)");
  if(areanorm==1){
    double dataarea = dataMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY());
    dataMnv = new MnvH2D(dataMnv->GetAreaNormalizedCopy());
    dataMnv->Scale(1e39, "width");
    mcMnv->Scale(1e39*dataarea/mcMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()), "width");   
    nommcMnv->Scale(1e39*dataarea/nommcMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    amuMnv->Scale(1e39*dataarea/amuMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    ncteqMnv->Scale(1e39*dataarea/ncteqMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    ncteqnuMnv->Scale(1e39*dataarea/ncteqnuMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    mkMnv->Scale(1e39*dataarea/mkMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    resminosMnv->Scale(1e39*dataarea/resminosMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    resnievesMnv->Scale(1e39*dataarea/resnievesMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    resminervaMnv->Scale(1e39*dataarea/resminervaMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    pion2p2hMnv->Scale(1e39*dataarea/pion2p2hMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    pionrpaMnv->Scale(1e39*dataarea/pionrpaMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    piontuneMnv->Scale(1e39*dataarea/piontuneMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    gibuuMnv->Scale(1e39*dataarea/gibuuMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    gibuu_v2021Mnv->Scale(1e39*dataarea/gibuu_v2021Mnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    nuwrosfMnv->Scale(1e39*dataarea/nuwrosfMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    nuwrolfgMnv->Scale(1e39*dataarea/nuwrolfgMnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    neutlfg105Mnv->Scale(1e39*dataarea/neutlfg105Mnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    neutsf103Mnv->Scale(1e39*dataarea/neutsf103Mnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
    //    neutsf121Mnv->Scale(1e39*dataarea/neutsf121Mnv->Integral(1,dataMnv->GetNbinsX(),1,dataMnv->GetNbinsY()),"width");
  }
  else{
    dataMnv->Scale(1e39, "width");
    mcMnv->Scale(1e39, "width");
    
    nommcMnv->Scale(1e39,"width");
    amuMnv->Scale(1e39,"width");
    ncteqMnv->Scale(1e39,"width");
    ncteqnuMnv->Scale(1e39,"width");
    mkMnv->Scale(1e39,"width");
    resminosMnv->Scale(1e39,"width");
    resnievesMnv->Scale(1e39,"width");
    resminervaMnv->Scale(1e39,"width");
    pion2p2hMnv->Scale(1e39,"width");
    pionrpaMnv->Scale(1e39,"width");
    piontuneMnv->Scale(1e39,"width");
    gibuuMnv->Scale(1e39,"width");
    gibuu_v2021Mnv->Scale(1e39,"width");
    nuwrosfMnv->Scale(1e39,"width");
    nuwrolfgMnv->Scale(1e39,"width");
    neutlfg105Mnv->Scale(1e39,"width");
    neutsf103Mnv->Scale(1e39,"width");
    //    neutsf121Mnv->Scale(1e39,"width");
  }


  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  

  TH2* nommc = new TH2D(nommcMnv->GetCVHistoWithStatError());
  TH2* amu= new TH2D(amuMnv->GetCVHistoWithStatError());
  TH2* ncteq= new TH2D(ncteqMnv->GetCVHistoWithStatError());
  TH2* ncteqnu= new TH2D(ncteqnuMnv->GetCVHistoWithStatError());
  TH2* mk= new TH2D(mkMnv->GetCVHistoWithStatError());
  TH2* resminos= new TH2D(resminosMnv->GetCVHistoWithStatError());
  TH2* resnieves= new TH2D(resnievesMnv->GetCVHistoWithStatError());
  TH2* resminerva= new TH2D(resminervaMnv->GetCVHistoWithStatError());
  TH2* pion2p2h= new TH2D(pion2p2hMnv->GetCVHistoWithStatError());
  TH2* pionrpa= new TH2D(pionrpaMnv->GetCVHistoWithStatError());
  TH2* piontune= new TH2D(piontuneMnv->GetCVHistoWithStatError());
  TH2* gibuu= new TH2D(gibuuMnv->GetCVHistoWithStatError());
  TH2* gibuu_v2021= new TH2D(gibuu_v2021Mnv->GetCVHistoWithStatError());
  TH2* nuwrosf= new TH2D(nuwrosfMnv->GetCVHistoWithStatError());
  TH2* nuwrolfg= new TH2D(nuwrolfgMnv->GetCVHistoWithStatError());
  TH2* neutlfg105 = new TH2D(neutlfg105Mnv->GetCVHistoWithStatError());
  TH2* neutsf103 = new TH2D(neutsf103Mnv->GetCVHistoWithStatError());
  //  TH2* neutsf121 = new TH2D(neutsf121Mnv->GetCVHistoWithStatError());




  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = MnvColors::GetColors(9);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);


  // These line and marker styles will be propagated to the 1D plots
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.5);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);

  dataStat->SetLineColor(kBlack);

  nommc->SetLineColor(mycolors[0]);  
  pion2p2h->SetLineColor(mycolors[1]);
  pionrpa->SetLineColor(mycolors[2]);
  piontune->SetLineColor(kViolet-3);

  amu->SetLineColor(mycolors[0]);
  ncteq->SetLineColor(mycolors[1]);
  ncteqnu->SetLineColor(mycolors[2]);

  mk->SetLineColor(mycolors[0]);
  resminos->SetLineColor(mycolors[1]);
  resnieves->SetLineColor(mycolors[2]);
  resminerva->SetLineColor(kViolet-3);

  gibuu->SetLineColor(mycolors[0]);//Green
  gibuu_v2021->SetLineColor(mycolors[1]);//
  nuwrosf->SetLineColor(mycolors[2]);
  nuwrolfg->SetLineColor(mycolors[6]);
  neutlfg105->SetLineColor(kViolet-3);
  neutsf103->SetLineColor(mycolors[7]);
  //  neutsf121->SetLineColor(mycolors[7]);


  amu->SetLineWidth(1.5);
  ncteq->SetLineWidth(1.5);
  ncteqnu->SetLineWidth(1.5);
  mk->SetLineWidth(1.5);
  resminos->SetLineWidth(1.5);
  resnieves->SetLineWidth(1.5);
  resminerva->SetLineWidth(1.5);
  pion2p2h->SetLineWidth(1.5);
  pionrpa->SetLineWidth(1.5);
  piontune->SetLineWidth(1.5);
  gibuu->SetLineWidth(1.5);
  gibuu_v2021->SetLineWidth(1.5);
  nommc->SetLineWidth(1.5);
  nuwrosf->SetLineWidth(1.5);
  nuwrolfg->SetLineWidth(1.5);
  neutlfg105->SetLineWidth(1.5);
  neutsf103->SetLineWidth(1.5);
  //  neutsf121->SetLineWidth(1.5);



  if(doRatio){
    //    TH2 *tmpden = (TH2*)mc->Clone("tmpden");
    //    tmpden->Sumw2(false);
    data->Divide(mc);
    dataStat->Divide(mc);
    amu->Divide(mc);
    ncteq->Divide(mc);
    ncteqnu->Divide(mc);
    mk->Divide(mc);
    resminos->Divide(mc);
    resnieves->Divide(mc);
    resminerva->Divide(mc);
    pion2p2h->Divide(mc);
    pionrpa->Divide(mc);
    piontune->Divide(mc);
    gibuu->Divide(mc);
    gibuu_v2021->Divide(mc);
    nommc->Divide(mc);
    nuwrosf->Divide(mc);
    nuwrolfg->Divide(mc);
    neutlfg105->Divide(mc);
    neutsf103->Divide(mc);
    //    neutsf121->Divide(mc);
    mc->Divide(mc);
  }
    

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
  histAndOpts.push_back(std::make_pair(dataStat, "histpe1"));
  histAndOpts.push_back(std::make_pair(mc,       "graph0LX"));
  if(modset ==0){//piontune,2p2h,rpa
    histAndOpts.push_back(std::make_pair(nommc,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(pion2p2h,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(pionrpa,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(piontune,       "graph0LX"));
  }
  else if(modset==1){
    histAndOpts.push_back(std::make_pair(amu,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(ncteq,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(ncteqnu,       "graph0LX"));
  }
  else if(modset==2){
    histAndOpts.push_back(std::make_pair(mk,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(resminos,       "graph0LX"));
    //    histAndOpts.push_back(std::make_pair(resnieves,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(resminerva,       "graph0LX"));
  }
  else if(modset==3){
    histAndOpts.push_back(std::make_pair(gibuu,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(gibuu_v2021,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(nuwrosf,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(nuwrolfg,       "graph0LX"));
    histAndOpts.push_back(std::make_pair(neutlfg105,   "graph0LX"));
    histAndOpts.push_back(std::make_pair(neutsf103,   "graph0LX"));
    //    histAndOpts.push_back(std::make_pair(neutsf121,   "graph0LX"));
  }

  histAndOpts.push_back(std::make_pair(data,     "histpe1"));



  // ----------------------------------------------------------------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers = GetScales(histAndOpts, true, 5,0.75);

  GridCanvas* gc=plotXAxis1DRebinPz(histAndOpts, "Muon Longitudinal Momentum (GeV/c)", "p_{t}", 4,4,800,500,doMultipliers ? &multipliers[0] : NULL);

  // Set the y range manually. Can also use gc->Remax() to guess automatically
  if(doRatio) gc->SetYLimits(0,1.99);
  else  gc->SetYLimits(0, 4.59);
  if(doRatio) gc->SetYTitle("Ratio to MINERvA Tune v1");
  else gc->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/(GeV/c)^{2}/c^{2}/Nucleon)");
  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.55, 0.05, 0.97, 0.32);
  leg->SetNColumns(2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MINERvA Tune v1", "l");

  if(modset ==0){//piontune,2p2h,rpa
    leg->AddEntry(nommc,"GENIE 2.12.6","l");
    leg->AddEntry(pion2p2h,"Low Recoil Enhancement","l");
    leg->AddEntry(pionrpa,"QE RPA","l");
    leg->AddEntry(piontune,"NonResPionTune Only","l");
  }
  else if(modset==1){
    leg->AddEntry(amu,"AMU DIS","l");
    leg->AddEntry(ncteq,"nCTEQ15 DIS","l");
    leg->AddEntry(ncteqnu,"nCTEQ#nu DIS","l");
  }
  else if(modset==2){
    leg->AddEntry(mk,"MK Model","l");
    leg->AddEntry(resminos,"Pion LowQ2 - MINOS","l");
    //    leg->AddEntry(resnieves,"Pion LowQ2 - Valencia","l");
    //    leg->AddEntry(resminerva,"Pion LowQ2 - MINERvA","l");
    leg->AddEntry(resminerva,"MINERvA Tune v2","l");
  }
  else if(modset==3){
    leg->AddEntry(gibuu,"GiBUU","l");
    leg->AddEntry(gibuu_v2021,"GiBUU v2021","l");
    leg->AddEntry(nuwrosf,"NuWro SF","l");
    leg->AddEntry(nuwrolfg,"NuWro LFG","l");
    leg->AddEntry(neutsf103,"NEUT SF","l");
    leg->AddEntry(neutlfg105,"NEUT LFG","l");
    //    leg->AddEntry(neutsf121,"NEUT SF M_{a}=1.21","l");
  }

  leg->Draw("SAME");

  if(doRatio){
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-pt-multiplier_ratio_models_set_%d_areanorm_%d.eps",modset,areanorm) : Form("nu-2d-xsec-comps-pt_ratio_models_set_%d_areanorm_%d.eps",modset,areanorm));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-pt-multiplier_ratio_models_set_%d_areanorm_%d.png",modset,areanorm) : Form("nu-2d-xsec-comps-pt_ratio_models_set_%d_areanorm_%d.png",modset,areanorm));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-pt-multiplier_ratio_models_set_%d_areanorm_%d.C",modset,areanorm) : Form("nu-2d-xsec-comps-pt_ratio_models_set_%d_areanorm_%d.C",modset,areanorm));
  }
  else{
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-pt-multiplier_models_set_%d_areanorm_%d.eps",modset,areanorm) : Form("nu-2d-xsec-comps-pt_models_set_%d_areanorm_%d.eps",modset,areanorm));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-pt-multiplier_models_set_%d_areanorm_%d.png",modset,areanorm) : Form("nu-2d-xsec-comps-pt_models_set_%d_areanorm_%d.png",modset,areanorm));
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-pt-multiplier_models_set_%d_areanorm_%d.C",modset,areanorm) : Form("nu-2d-xsec-comps-pt_models_set_%d_areanorm_%d.C",modset,areanorm));
  }
  

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
  vector<double> multipliers2 = GetScales(histAndOpts, false, 2.5,0.75);

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  //GridCanvas* gc2=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV/c)", "p_{||}",4,4,800,500,doMultipliers ? &multipliers2[0] : NULL);
  GridCanvas* gc2=plotYAxis1DRebinPt(histAndOpts, "Muon Transverse Momentum (GeV/c)", "p_{||}",4,4,800,500,doMultipliers ? &multipliers2[0] : NULL);
  if(doRatio) gc2->SetYLimits(0,1.99);
  else  gc2->SetYLimits(0, 2.49);
  if(doRatio) gc2->SetYTitle("Ratio to MINERvA Tune v1");
  else gc2->SetYTitle("d^{2}#sigma/dp_{T}dp_{||} (x10^{-39} cm^{2}/(GeV/c)^{2}/c^{2}/Nucleon)");
  gc2->Modified();
  if(doRatio){
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-pz-multiplier_ratio_models_set_%d_areanorm_%d.eps",modset,areanorm) : Form("nu-2d-xsec-comps-pz_ratio_models_set_%d_areanorm_%d.eps",modset,areanorm));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-pz-multiplier_ratio_models_set_%d_areanorm_%d.png",modset,areanorm) : Form("nu-2d-xsec-comps-pz_ratio_models_set_%d_areanorm_%d.png",modset,areanorm));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-pz-multiplier_ratio_models_set_%d_areanorm_%d.C",modset,areanorm) : Form("nu-2d-xsec-comps-pz_ratio_models_set_%d_areanorm_%d.C",modset,areanorm));
  }
  else{
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-pz-multiplier_models_set_%d_areanorm_%d.eps",modset,areanorm) : Form("nu-2d-xsec-comps-pz_models_set_%d_areanorm_%d.eps",modset,areanorm));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-pz-multiplier_models_set_%d_areanorm_%d.png",modset,areanorm) : Form("nu-2d-xsec-comps-pz_models_set_%d_areanorm_%d.png",modset,areanorm));
    gc2->Print(doMultipliers ? Form("nu-2d-xsec-comps-pz-multiplier_models_set_%d_areanorm_%d.C",modset,areanorm) : Form("nu-2d-xsec-comps-pz_models_set_%d_areanorm_%d.C",modset,areanorm));
  }

}

int main(int argc, char* argv[])
{
  //  makePlots(true,true,argv[1]); //multipliers on ratio don't make sense
  for(int i=0;i<4;i++){
    for(int j=0;j<2;j++){
      makePlots(true,false,argv[1],i,j);
      makePlots(false,true,argv[1],i,j);
      makePlots(false,false,argv[1],i,j);
    }
  }
  return 0;
}
