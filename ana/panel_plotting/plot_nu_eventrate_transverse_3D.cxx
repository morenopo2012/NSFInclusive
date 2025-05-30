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
#include "/minerva/app/users/drut1186/cmtuser/Minerva_v22r1p1_CCQENu/Ana/CCQENu/ana_common/src/CCQENuBinning.cxx"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;
using namespace CCQENU_ANA;

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


void makePlots(bool doMultipliers,bool doGenies,string location,string varx, string vary, bool do3D)
{
  CCQENuBinning *minmodbinner = new CCQENuBinning();
  //Binning needed for 3D plotting
  //  axis_binning dthetaPerpbins_TKI = neutbinner->GetThetaPerpBinsDegree_TKI();
  //  axis_binning dthetaReactbins_TKI = neutbinner->GetThetaReactionPlaneBinsDegree_TKI();
  axis_binning muonPtbins     = minmodbinner->GetMuonPtBinsGeV();//13
  axis_binning muonPzbins     = minmodbinner->GetMuonPzBinsGeV();//13
  axis_binning dalphatbins_TKI    = minmodbinner->GetDalphaT_TKI();
  axis_binning dphitbins_TKI    = minmodbinner->GetDphiT_TKI();
  axis_binning pnbins_TKI    = minmodbinner->GetPn_TKI();
  axis_binning dptbins_TKI    = minmodbinner->GetDeltaPt_TKI();
  axis_binning dptxbins_TKI    = minmodbinner->GetDeltaPtx_TKI();
  axis_binning dptybins_TKI    = minmodbinner->GetDeltaPty_TKI();
  axis_binning protonThetabins_TKI = minmodbinner->GetProtonThetaBins_TKI();
  axis_binning protonKineticbins_TKI = minmodbinner->GetProtonEnergyBinsGeV_TKI();
  axis_binning signeddalphatbins_TKI = minmodbinner->GetSignedDeltaAlphaT_TKI(); 
  axis_binning signeddphitbins_TKI = minmodbinner->GetSignedDeltaPhiT_TKI();
  axis_binning recoil_TKI = minmodbinner->Get3DRecoilBins();
  axis_binning signbins          = minmodbinner->GetSignedBins();


  map<int,axis_binning > axis_bin_map;
  map<string,int > axis_name;
  axis_bin_map[0] = muonPzbins;
  axis_bin_map[1] = dalphatbins_TKI;
  axis_bin_map[2] = dphitbins_TKI;
  axis_bin_map[3] = pnbins_TKI;
  axis_bin_map[4] = dptbins_TKI;
  axis_bin_map[5] = dptxbins_TKI;
  axis_bin_map[6] = dptybins_TKI;
  axis_bin_map[7] = protonKineticbins_TKI;
  axis_bin_map[8] = protonThetabins_TKI;
  axis_bin_map[9] = signbins;
  axis_bin_map[10] = signeddalphatbins_TKI;
  axis_bin_map[11] = signeddphitbins_TKI;
  //  axis_bin_map[12] = dthetaReactbins_TKI;
  //  axis_bin_map[13] = dthetaPerpbins_TKI;
  axis_bin_map[14] = dptybins_TKI;
  axis_bin_map[15] = dptxbins_TKI;
  axis_bin_map[16] = dptybins_TKI;
  axis_bin_map[17] = dptxbins_TKI;
  axis_bin_map[18] = recoil_TKI; 


  axis_name["muonPz"] = 0;
  axis_name["dalphat"] = 1;
  axis_name["dphit"] = 2;
  axis_name["pn"] = 3;
  axis_name["dpt"] = 4;
  axis_name["dptx"] = 5;
  axis_name["dpty"] = 6;
  axis_name["tp"] = 7;
  axis_name["ptheta"] = 8;
  axis_name["signed"] = 9;
  axis_name["signeddalphat"] = 10;
  axis_name["signeddphit"] = 11;
  axis_name["dthetaR"] = 12;
  axis_name["dthetaP"] = 13;
  axis_name["dPR"] = 14;
  axis_name["dPP"] = 15;
  axis_name["dPRi"] = 16;
  axis_name["dPPi"] = 17;
  axis_name["recoil"] = 18;

  int varx_int = axis_name[varx];
  int vary_int = axis_name[vary];

  vector<vector<double> >full3D;
  full3D.push_back(axis_bin_map[varx_int].bin_edges);
  full3D.push_back(muonPtbins.bin_edges);
  full3D.push_back(axis_bin_map[vary_int].bin_edges);

  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);  
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  string var = varx+"_"+vary;
  string varspecfilename = var;
  string varspec = var;
  if(vary=="ptmu") varspecfilename=varx;
  cout << varspecfilename << endl;
  string filename = Form("%s_CV_ISI_NuWro_LFG//selection_Signal/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_%s_CombinedPlaylists.root",location.c_str(),varspecfilename.c_str());//Final result
  TFile f1(filename.c_str());

  MnvH2D* dataMnv, *mcMnv, *mcMnv_qelike_qe, *mcMnv_qelike_res, *mcMnv_qelike_dis, *mcMnv_qelike_2p2h, *mcMnv_qelikenot;

  if(do3D){
    dataMnv = (MnvH2D*)f1.Get(Form("h_%s_data",varspec.c_str()));
    mcMnv = (MnvH2D*)f1.Get(Form("h_%s_mc",varspec.c_str()));
    
    mcMnv_qelike_qe = (MnvH2D*)f1.Get(Form("h_%s_qelike_qe",varspec.c_str()));//Get from N track
    mcMnv_qelike_res = (MnvH2D*)f1.Get(Form("h_%s_qelike_res",varspec.c_str()));//Get from N track
    mcMnv_qelike_dis = (MnvH2D*)f1.Get(Form("h_%s_qelike_dis",varspec.c_str()));//Get from N track
    mcMnv_qelike_2p2h = (MnvH2D*)f1.Get(Form("h_%s_qelike_2p2h",varspec.c_str()));//Get from N track
    mcMnv_qelikenot   = (MnvH2D*)f1.Get(Form("h_%s_qelikenot",varspec.c_str()));//Get from N track
  
  }
  else{
    dataMnv=(MnvH2D*)f1.Get(Form("h_%s_data_remapped",varspec.c_str()));
    mcMnv=(MnvH2D*)f1.Get(Form("h_%s_mc_remapped",varspec.c_str()));
    
    mcMnv_qelike_qe = (MnvH2D*)f1.Get(Form("h_%s_qelike_qe_remapped",varspec.c_str()));//Get from N track
    mcMnv_qelike_res = (MnvH2D*)f1.Get(Form("h_%s_qelike_res_remapped",varspec.c_str()));//Get from N track
    mcMnv_qelike_dis = (MnvH2D*)f1.Get(Form("h_%s_qelike_dis_remapped",varspec.c_str()));//Get from N track
    mcMnv_qelike_2p2h = (MnvH2D*)f1.Get(Form("h_%s_qelike_2p2h_remapped",varspec.c_str()));//Get from N track
    mcMnv_qelikenot   = (MnvH2D*)f1.Get(Form("h_%s_qelikenot_remapped",varspec.c_str()));//Get from N track
  }

  vector<TH2D*> data,dataStat, mc, mc_qelike_qe, mc_qelike_res, mc_qelike_2p2h, mc_qelike_dis, mc_qelikenot;

  
  if(do3D){
    cout << "Doing data project" << endl;
    data = my3d->Get2DHistos(dataMnv,true);
    dataStat = my3d->Get2DHistos(dataMnv,false);
    cout << "Doing mc project" << endl;
    mc =my3d->Get2DHistos(mcMnv,false);
    cout << "Doing mc qelike project" << endl;
    mc_qelike_qe =my3d->Get2DHistos(mcMnv_qelike_qe,false);
    mc_qelike_res =my3d->Get2DHistos(mcMnv_qelike_res,false);
    mc_qelike_2p2h =my3d->Get2DHistos(mcMnv_qelike_2p2h,false);
    mc_qelike_dis =my3d->Get2DHistos(mcMnv_qelike_dis,false);
    cout << "Doing mc_qelikenot project" << endl;
    mc_qelikenot =my3d->Get2DHistos(mcMnv_qelikenot,false);
  }
  else{
    // Get the data histogram with stat error and with total error
    // separately so we can plot them both for inner and outer ticks
    dataStat.push_back(new TH2D(dataMnv->GetCVHistoWithStatError()));
    data.push_back(new TH2D(dataMnv->GetCVHistoWithError()));
    mc.push_back(new TH2D(mcMnv->GetCVHistoWithStatError()));
    
    mc_qelike_qe.push_back(new TH2D(mcMnv_qelike_qe->GetCVHistoWithStatError()));
    mc_qelike_res.push_back(new TH2D(mcMnv_qelike_res->GetCVHistoWithStatError()));
    mc_qelike_dis.push_back(new TH2D(mcMnv_qelike_dis->GetCVHistoWithStatError()));
    mc_qelike_2p2h.push_back(new TH2D(mcMnv_qelike_2p2h->GetCVHistoWithStatError()));
    mc_qelikenot.push_back(new TH2D(mcMnv_qelikenot->GetCVHistoWithStatError()));
  }



  for(int i=0;i<data.size();i++){
    cout << "Doing " << i << endl;
    if(do3D){
      double zwidth = axis_bin_map[vary_int].bin_edges[i]-axis_bin_map[vary_int].bin_edges[i-1];
      
      
      //Move this later
      data[i]->Scale(1-6*zwidth, "width");
      dataStat[i]->Scale(1-6*zwidth, "width");
      mc[i]->Scale(1-6*zwidth, "width");
      mc_qelike_qe[i]->Scale(1-6*zwidth,"width");
      mc_qelike_res[i]->Scale(1-6*zwidth,"width");
      mc_qelike_dis[i]->Scale(1-6*zwidth,"width");
      mc_qelike_2p2h[i]->Scale(1-6*zwidth,"width");
      mc_qelikenot[i]->Scale(1-6*zwidth,"width");
    }
    else{
      data[i]->Scale(1-6, "width");
      dataStat[i]->Scale(1-6, "width");
      mc[i]->Scale(1-6, "width");
      mc_qelike_qe[i]->Scale(1-6,"width");
      mc_qelike_res[i]->Scale(1-6,"width");
      mc_qelike_dis[i]->Scale(1-6,"width");
      mc_qelike_2p2h[i]->Scale(1-6,"width");
      mc_qelikenot[i]->Scale(1-6,"width");
    }

    
    
    
    
    
    
    // These line and marker styles will be propagated to the 1D plots
    vector<int> mycolors = getColors(2);
    mc[i]->SetLineColor(kRed);
    mc[i]->SetLineWidth(2);
    mc_qelike_qe[i]->SetLineColor(mycolors[3]);
    mc_qelike_res[i]->SetLineColor(mycolors[4]);
    mc_qelike_dis[i]->SetLineColor(mycolors[5]);
    mc_qelike_2p2h[i]->SetLineColor(mycolors[16]);
    mc_qelikenot[i]->SetLineColor(mycolors[8]);

    // These line and marker styles will be propagated to the 1D plots
    data[i]->SetMarkerStyle(kFullCircle);
    data[i]->SetMarkerSize(0.5);
    data[i]->SetLineColor(kBlack);
    data[i]->SetLineWidth(2);

    dataStat[i]->SetLineColor(kBlack);


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
    histAndOpts.push_back(std::make_pair(dataStat[i], "histpe1"));
    histAndOpts.push_back(std::make_pair(mc[i],       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_qe[i],       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res[i],       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_dis[i],       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h[i],       "hist l"));
    histAndOpts.push_back(std::make_pair(mc_qelikenot[i],       "hist l"));
    histAndOpts.push_back(std::make_pair(data[i],     "histpe1"));



    // ----------------------------------------------------------------------------------
    //
    // First make pt in bins of pz

    // Values to multiply each bin by to get them on a similar range
    vector<double> multi_1 = GetScales(histAndOpts, false);

    GridCanvas* gc=NULL;
    if(do3D) gc= plotYAxis1D(histAndOpts, vary.c_str() ,varx.c_str() ,doMultipliers ? &multi_1[0] : NULL);
    else  gc= plotYAxis1D(histAndOpts, "Muon Transverse Momentum" ,varx.c_str() ,doMultipliers ? &multi_1[0] : NULL);
  
    string mystring =     Form("%.2f < %s [GeV] < %.2f",axis_bin_map[vary_int].bin_edges[i-1],vary.c_str(),axis_bin_map[vary_int].bin_edges[i]);
    TLatex *mytex = new TLatex();
    mytex->SetTextSize(0.05);
    mytex->SetTextFont(42);

    // Set the y range manually. Can also use gc->Remax() to guess automatically
    gc->SetYLimits(0, 4.99);
  

    //Label thy axis
    gc->SetYTitle("Event Rate(*1-6) per GeV^{2}");
    gc->SetXTitle(vary.c_str());
    if(do3D)gc->SetXTitle("Muon Transverse");
    
    gc->Modified();
    // Example of adding a legend. The co-ordinate system is NDC on the
    // entire canvas, ie (0,0) in the bottom left corner of the canvas
    // (not the individual pad), and (1,1) in the top right
    //  TLegend* leg=new TLegend(0.17, 0.7, 0.31, 0.9);
    TLegend* leg=new TLegend(0.4, 0.1, 0.9, 0.3);
    leg->SetNColumns(2);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(data[i], "MINERvA data", "lpe");
    leg->AddEntry(mc[i], "MnvGENIEv1.0.1", "l");
    leg->AddEntry(mc_qelike_qe[i],"QE","l");
    leg->AddEntry(mc_qelike_res[i],"Resonant","l");
    leg->AddEntry(mc_qelike_dis[i],"DIS","l");
    leg->AddEntry(mc_qelike_2p2h[i],"2p2h","l");
    leg->AddEntry(mc_qelikenot[i],"Background","l");
    leg->Draw("SAME");
    if(do3D)mytex->DrawLatex(0.35,0.96,mystring.c_str());
    if(do3D){
      gc->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projy-%s-bin-%d-multiplier.eps",var.c_str(),i) : Form("nu-2d-evtrate-comps-projy-%s-bin-%d.eps",var.c_str(),i));
      gc->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projy-%s-bin-%d-multiplier.png",var.c_str(),i) : Form("nu-2d-evtrate-comps-projy-%s-bin-%d.png",var.c_str(),i));
      gc->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projy-%s-bin-%d-multiplier.C",var.c_str(),i) : Form("nu-2d-evtrate-comps-projy-%s-bin-%d.C",var.c_str(),i));
    }
    else{
      gc->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projy-%s-multiplier.eps",var.c_str()) : Form("nu-2d-evtrate-comps-projy-%s.eps",var.c_str()));
      gc->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projy-%s-multiplier.png",var.c_str()) : Form("nu-2d-evtrate-comps-projy-%s.png",var.c_str()));
      gc->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projy-%s-multiplier.C",var.c_str()) : Form("nu-2d-evtrate-comps-projy-%s.C",var.c_str()));
    }


    // ----------------------------------------------------------------------------------
    //
    // Now make pz in bins of pt. It's all the same

    vector<double> multi_2 = GetScales(histAndOpts,true);
  


    // Values to multiply each bin by to get them on a similar range


    // plotpz1D fiddles the x axis values to squash up the tail so it
    // doesn't take up all the horizontal space.
    GridCanvas* gc2=NULL;
    if(do3D) gc2 =plotXAxis1D(histAndOpts, varx.c_str() , vary.c_str() ,doMultipliers ? &multi_2[0] : NULL);
    else  gc2 =plotXAxis1D(histAndOpts, varx.c_str() , vary.c_str() ,doMultipliers ? &multi_2[0] : NULL);
    gc2->SetYLimits(0, 4.99);
    gc2->SetYTitle("Event rate(x1-6) per GeV^{2}");
    gc2->SetXTitle(varx.c_str());
    if(do3D)    gc2->SetXTitle(varx.c_str());
    leg->Draw("SAME");
    if(do3D)mytex->DrawLatex(0.35,0.96,mystring.c_str());
    gc2->Modified();
    if(do3D){
      gc2->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projx-%s-bin-%d-multiplier.eps",var.c_str(),i) : Form("nu-2d-evtrate-comps-projx-%s-bin-%d.eps",var.c_str(),i));
      gc2->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projx-%s-bin-%d-multiplier.png",var.c_str(),i) : Form("nu-2d-evtrate-comps-projx-%s-bin-%d.png",var.c_str(),i));
      gc2->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projx-%s-bin-%d-multiplier.C",var.c_str(),i) : Form("nu-2d-evtrate-comps-projx-%s-bin-%d.C",var.c_str(),i));
    }
    else{
      gc2->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projx-%s-multiplier.eps",var.c_str()) : Form("nu-2d-evtrate-comps-projx-%s.eps",var.c_str()));
      gc2->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projx-%s-multiplier.png",var.c_str()) : Form("nu-2d-evtrate-comps-projx-%s.png",var.c_str()));
      gc2->Print(doMultipliers ? Form("nu-2d-evtrate-comps-projx-%s-multiplier.C",var.c_str()) : Form("nu-2d-evtrate-comps-projx-%s.C",var.c_str()));
    }
  }
}

int main(int argc, char* argv[])
{
  // makePlots(true,true,argv[1],argv[2]);
  // makePlots(true,false,argv[1],argv[2]);
  // makePlots(false,true,argv[1],argv[2]);
  //  makePlots(false,false,argv[1],argv[2],argv[3]);
  //makePlots(true,false,argv[1],argv[2],argv[3],false);
    makePlots(true,false,argv[1],argv[2],argv[3],true);

  return 0;
}
