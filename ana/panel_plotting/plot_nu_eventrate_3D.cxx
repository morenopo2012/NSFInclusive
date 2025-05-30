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

void makePlots(bool doMultipliers, bool doTracks, bool doBkgSub, string location, bool withlowrecoilremoved, bool withresfsi, bool withqefsi, bool doRecoilFit, bool doRatio, bool doZoom,bool domcerr,bool doconst)
{

  cout << "Running with " << withlowrecoilremoved << "\t" << withresfsi << "\t" << withqefsi  << "\t" << domcerr << "\t" << doconst<< endl;
  //  return;
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
  pt3Dbins.push_back(0.075);//added ME
  pt3Dbins.push_back(0.15);
  pt3Dbins.push_back(0.25);//added ME
  pt3Dbins.push_back(0.325);//added ME
  pt3Dbins.push_back(0.4);
  pt3Dbins.push_back(0.475);//added ME
  pt3Dbins.push_back(0.55);//added ME
  pt3Dbins.push_back(0.7);//added ME
  pt3Dbins.push_back(0.85);//added ME
  pt3Dbins.push_back(1.0);

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


  /*
  recoil3Dbins.push_back(0);
  recoil3Dbins.push_back(0.02);
  for(int i=0;i<9;i++)recoil3Dbins.push_back(i*0.04+0.04);
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*0.200+0.400);
  recoil3Dbins.push_back(0.7999); 
    
  */

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




  //three files 1-track, 2+track, N-track
  //CV
  TFile f1(Form("%s_CV_Zexp_Bubble/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_CombinedPlaylists.root",location.c_str()));//1track
  TFile f2(Form("%s_CV_Zexp_Bubble/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_CombinedPlaylists.root",location.c_str()));//2track
  TFile f3(Form("%s_CV_Zexp_Bubble/MuonEventSelection_MakeFlux-1_Multiplicity-0_Sample-Signal_CombinedPlaylists.root",location.c_str()));//Ntrack
  //CV without low recoil tune
  TFile f4(Form("%s_pion_rpa/MuonEventSelection_MakeFlux-1_Multiplicity-0_Sample-Signal_CombinedPlaylists.root",location.c_str()));//NTrack

  //need pzmuptmu
  MnvH2D* mcMnv=(MnvH2D*)f3.Get("h_pzptrec_mc");//Get from N track
  MnvH2D* dataMnv = (MnvH2D*)f3.Get("h_pzptrec_data");//Get from N track
  MnvH2D* SingleTrackMnv = (MnvH2D*)f1.Get("h_pzptrec_data");//Get from 1 track
  MnvH2D* MultiTrackMnv = (MnvH2D*)f2.Get("h_pzptrec_data");//Get from 2+ track
  MnvH2D* SingleTrackMCMnv = (MnvH2D*)f1.Get("h_pzptrec_mc");//Get from 1 track
  MnvH2D* MultiTrackMCMnv = (MnvH2D*)f2.Get("h_pzptrec_mc");//Get from 2+ track

  MnvH2D* mcMnv_qelike_qe = (MnvH2D*)f3.Get("h_pzptrec_qelike_qe");//Get from N track
  MnvH2D* mcMnv_qelike_res = (MnvH2D*)f3.Get("h_pzptrec_qelike_res");//Get from N track
  MnvH2D* mcMnv_qelike_dis = (MnvH2D*)f3.Get("h_pzptrec_qelike_dis");//Get from N track
  MnvH2D* mcMnv_qelike_2p2h = (MnvH2D*)f3.Get("h_pzptrec_qelike_2p2h");//Get from N track
  MnvH2D* mcMnv_qelike_qe_proton_fsi = (MnvH2D*)f3.Get("h_pzptrec_qelike_qe_proton_fsi");//Get from N track
  MnvH2D* mcMnv_qelike_qe_neutron_fsi = (MnvH2D*)f3.Get("h_pzptrec_qelike_qe_neutron_fsi");//Get from N track
  MnvH2D* mcMnv_qelike_res_proton_fsi = (MnvH2D*)f3.Get("h_pzptrec_qelike_res_proton_fsi");//Get from N track
  MnvH2D* mcMnv_qelike_res_neutron_fsi = (MnvH2D*)f3.Get("h_pzptrec_qelike_res_neutron_fsi");//Get from N track

  MnvH2D* mcMnv_qelike_2p2h_nntune = new MnvH2D(*mcMnv->GetVertErrorBand("Low_Recoil_2p2h_Tune")->GetHist(0));
  MnvH2D* mcMnv_qelike_2p2h_nptune = new MnvH2D(*mcMnv->GetVertErrorBand("Low_Recoil_2p2h_Tune")->GetHist(1));
  MnvH2D* mcMnv_qelike_2p2h_qetune = new MnvH2D(*mcMnv->GetVertErrorBand("Low_Recoil_2p2h_Tune")->GetHist(2));


  MnvH2D* mcMnv_qelike_2p2h_no_lowrec = (MnvH2D*)f4.Get("h_pzptrec_qelike_2p2h");//Get from N track


  //bkgs
  MnvH2D* mcMnv_1track_bkg = (MnvH2D*)f1.Get("h_pzptrec_qelikenot");//Get from 1 track
  MnvH2D* mcMnv_2track_bkg = (MnvH2D*)f2.Get("h_pzptrec_qelikenot");//Get from 2+ track
  
  TFile f_bkgscale(Form("%s_CV_Zexp_Bubble/BkgWeights3D.root",location.c_str()));
  MnvH2D *scalehist = (MnvH2D*)f_bkgscale.Get("bkgconstants");


  MnvH2D* mcMnv_bkg = (MnvH2D*)mcMnv_1track_bkg->Clone("tmp");
  mcMnv_bkg->Add(mcMnv_2track_bkg);

  MnvH2D* mcMnv_bkg_constrained = (MnvH2D*)mcMnv_bkg->Clone("const_tmp");
  mcMnv_bkg_constrained->Multiply(mcMnv_bkg_constrained,scalehist);


  if(doconst){
    mcMnv->Add(mcMnv_bkg,-1.0);//subtract off original bkg
    mcMnv->Add(mcMnv_bkg_constrained,1.0);//add new constrained on
    mcMnv_bkg = (MnvH2D*)mcMnv_bkg_constrained->Clone("constrained_bkg");//Make the bkg plotted the constrained...
  }
  
  

  
  dataMnv->GetXaxis()->SetTitle("p_{||} (GeV/c)");
  /*
  if(doBkgSub){
    dataMnv->Add(mcMnv_bkg,-1);
    mcMnv->Add(mcMnv_bkg,-1);
    SingleTrackMnv->Add(mcMnv_1track_bkg,-1);
    MultiTrackMnv->Add(mcMnv_2track_bkg,-1);
    SingleTrackMCMnv->Add(mcMnv_1track_bkg,-1);
    MultiTrackMCMnv->Add(mcMnv_2track_bkg,-1);
  }
  */

  std::vector<TH2D*> dataStat = my3d->Get2DHistos(dataMnv,false);
  std::vector<TH2D*> data = my3d->Get2DHistos(dataMnv,true);
  std::vector<TH2D*> mc = my3d->Get2DHistos(mcMnv,false);
  std::vector<TH2D*> mc_witherr = my3d->Get2DHistos(mcMnv,true);
  std::vector<TH2D*> mc_qelike_qe = my3d->Get2DHistos(mcMnv_qelike_qe,false);
  std::vector<TH2D*> mc_qelike_res = my3d->Get2DHistos(mcMnv_qelike_res,false);
  std::vector<TH2D*> mc_qelike_dis = my3d->Get2DHistos(mcMnv_qelike_dis,false);
  std::vector<TH2D*> mc_qelike_2p2h = my3d->Get2DHistos(mcMnv_qelike_2p2h,false);
  std::vector<TH2D*> mc_qelike_2p2h_no_lowrec = my3d->Get2DHistos(mcMnv_qelike_2p2h_no_lowrec,false);
  std::vector<TH2D*> mc_qelike_qe_proton_fsi = my3d->Get2DHistos(mcMnv_qelike_qe_proton_fsi,false);
  std::vector<TH2D*> mc_qelike_qe_neutron_fsi = my3d->Get2DHistos(mcMnv_qelike_qe_neutron_fsi,false);
  std::vector<TH2D*> mc_qelike_res_proton_fsi = my3d->Get2DHistos(mcMnv_qelike_res_proton_fsi,false);
  std::vector<TH2D*> mc_qelike_res_neutron_fsi = my3d->Get2DHistos(mcMnv_qelike_res_neutron_fsi,false);

  std::vector<TH2D*> mc_qelike_2p2h_nptune = my3d->Get2DHistos(mcMnv_qelike_2p2h_nptune,false);
  std::vector<TH2D*> mc_qelike_2p2h_nntune = my3d->Get2DHistos(mcMnv_qelike_2p2h_nntune,false);
  std::vector<TH2D*> mc_qelike_2p2h_qetune = my3d->Get2DHistos(mcMnv_qelike_2p2h_qetune,false);


  std::vector<TH2D*> mc_bkg = my3d->Get2DHistos(mcMnv_bkg,false);
  /*
  std::vector<TH2D*> data_1track = my3d->Get2DHistos(SingleTrackMnv,false);
  std::vector<TH2D*> data_2track = my3d->Get2DHistos(MultiTrackMnv,false);
  std::vector<TH2D*> mc_1track = my3d->Get2DHistos(SingleTrackMCMnv,false);
  std::vector<TH2D*> mc_2track = my3d->Get2DHistos(MultiTrackMCMnv,false);
  */

  // // Get the data histogram with stat error and with total error
  // // separately so we can plot them both for inner and outer ticks
  // TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  // TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  // TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  // TH2* mc_qelike_qe = new TH2D(mcMnv_qelike_qe->GetCVHistoWithStatError());
  // TH2* mc_qelike_res = new TH2D(mcMnv_qelike_res->GetCVHistoWithStatError());
  // TH2* mc_qelike_dis = new TH2D(mcMnv_qelike_dis->GetCVHistoWithStatError());
  // TH2* mc_qelike_2p2h = new TH2D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());
  // //  TH2* mc_qelike_2p2h_no_lowrec = new TH2D(mcMnv_qelike_2p2h_no_lowrec->GetCVHistoWithStatError());
  // TH2* mc_bkg = new TH2D(mcMnv_bkg->GetCVHistoWithStatError());
  // TH2* data_1track = new TH2D(SingleTrackMnv->GetCVHistoWithStatError());
  // TH2* data_2track = new TH2D(MultiTrackMnv->GetCVHistoWithStatError());
  // TH2* mc_1track = new TH2D(SingleTrackMCMnv->GetCVHistoWithStatError());
  // TH2* mc_2track = new TH2D(MultiTrackMCMnv->GetCVHistoWithStatError());
  for(int i=1;i<pz3Dbins.size();i++){

    double pzwidth = pz3Dbins[i]-pz3Dbins[i-1];

    dataStat[i]->Scale(1e-5/pzwidth,"width");
    data[i]->Scale(1e-5/pzwidth,"width");
    mc[i]->Scale(1e-5/pzwidth,"width");
    mc_witherr[i]->Scale(1e-5/pzwidth,"width");
    mc_qelike_qe[i]->Scale(1e-5/pzwidth,"width");
    mc_qelike_res[i]->Scale(1e-5/pzwidth,"width");
    mc_qelike_dis[i]->Scale(1e-5/pzwidth,"width");
    mc_qelike_2p2h[i]->Scale(1e-5/pzwidth,"width");
    mc_qelike_qe_proton_fsi[i]->Scale(1e-5/pzwidth,"width");
    mc_qelike_qe_neutron_fsi[i]->Scale(1e-5/pzwidth,"width");
    mc_qelike_res_proton_fsi[i]->Scale(1e-5/pzwidth,"width");
    mc_qelike_res_neutron_fsi[i]->Scale(1e-5/pzwidth,"width");

    mc_qelike_2p2h_nptune[i]->Scale(1e-5/pzwidth,"width");
    mc_qelike_2p2h_nntune[i]->Scale(1e-5/pzwidth,"width");
    mc_qelike_2p2h_qetune[i]->Scale(1e-5/pzwidth,"width");


    mc_bkg[i]->Scale(1e-5/pzwidth,"width");
    /*
    data_1track[i]->Scale(1e-5/pzwidth,"width");
    data_2track[i]->Scale(1e-5/pzwidth,"width");
    mc_1track[i]->Scale(1e-5/pzwidth,"width");
    mc_2track[i]->Scale(1e-5/pzwidth,"width");
    */
    mc_qelike_2p2h_no_lowrec[i]->Scale(1e-5/pzwidth,"width");
    

    if(doRatio){

    dataStat[i]->Divide(mc[i]);
    data[i]->Divide(mc[i]);
    mc_qelike_qe[i]->Divide(mc[i]);
    mc_qelike_res[i]->Divide(mc[i]);
    mc_qelike_dis[i]->Divide(mc[i]);
    mc_qelike_2p2h[i]->Divide(mc[i]);
    mc_qelike_qe_proton_fsi[i]->Divide(mc[i]);
    mc_qelike_qe_neutron_fsi[i]->Divide(mc[i]);
    mc_qelike_res_proton_fsi[i]->Divide(mc[i]);
    mc_qelike_res_neutron_fsi[i]->Divide(mc[i]);
    mc_qelike_2p2h_nptune[i]->Divide(mc[i]);
    mc_qelike_2p2h_nntune[i]->Divide(mc[i]);
    mc_qelike_2p2h_qetune[i]->Divide(mc[i]);

    mc_bkg[i]->Divide(mc[i]);
    /*
    data_1track[i]->Divide(mc[i]);
    data_2track[i]->Divide(mc[i]);
    mc_1track[i]->Divide(mc[i]);
    mc_2track[i]->Divide(mc[i]);
    */
    mc_qelike_2p2h_no_lowrec[i]->Divide(mc[i]);
    mc_witherr[i]->Divide(mc[i]);
    mc[i]->Divide(mc[i]);

    }
    

    // These line and marker styles will be propagated to the 1D plots
    //vector<int> mycolors = getColors(2);
    //    vector<int> mycolors2 = getColors(1);
    vector<int> mycolors = MnvColors::GetColors(9);
    vector<int> mycolors2 = MnvColors::GetColors(6);

    mc[i]->SetLineColor(kRed);
    //    if(domcerr) mc[i]->SetLineColor(kRed+2);
    mc[i]->SetLineWidth(2);

    mc_witherr[i]->SetLineColor(kRed);
    mc_witherr[i]->SetLineWidth(2);
    mc_witherr[i]->SetFillColor(kRed);
    mc_witherr[i]->SetFillStyle(3002);
    /*    
    mc_1track[i]->SetLineColor(kViolet);
    mc_1track[i]->SetLineStyle(2);
    mc_2track[i]->SetLineColor(kViolet);
    mc_2track[i]->SetLineStyle(1);
    */
    //need to add signal and bkg colors
    mc_qelike_qe[i]->SetLineColor(mycolors[2]);
    mc_qelike_qe_proton_fsi[i]->SetLineColor(mycolors[2]);
    mc_qelike_qe_neutron_fsi[i]->SetLineColor(mycolors[2]);
    mc_qelike_qe_proton_fsi[i]->SetLineStyle(2);
    mc_qelike_qe_neutron_fsi[i]->SetLineStyle(9);
    mc_qelike_res[i]->SetLineColor(mycolors[1]);
    mc_qelike_res_proton_fsi[i]->SetLineColor(mycolors[1]);
    mc_qelike_res_neutron_fsi[i]->SetLineColor(mycolors[1]);
    mc_qelike_res_proton_fsi[i]->SetLineStyle(2);
    mc_qelike_res_neutron_fsi[i]->SetLineStyle(9);

    mc_qelike_2p2h_nptune[i]->SetLineColor(mycolors[1]);
    mc_qelike_2p2h_nntune[i]->SetLineColor(mycolors[2]);
    mc_qelike_2p2h_qetune[i]->SetLineColor(mycolors2[5]);

    mc_qelike_dis[i]->SetLineColor(mycolors[7]);
    mc_qelike_2p2h[i]->SetLineColor(mycolors[0]);
    mc_qelike_2p2h_no_lowrec[i]->SetLineColor(mycolors[0]);
    mc_qelike_2p2h_no_lowrec[i]->SetLineStyle(2);
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
    /*
    data_1track[i]->SetLineColor(kBlack);
    data_1track[i]->SetLineStyle(2);
    data_2track[i]->SetLineColor(kBlack);
    data_2track[i]->SetLineStyle(1);
    */
    
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
    histAndOpts.push_back(std::make_pair(data[i],     "histpe"));
    if(!domcerr) histAndOpts.push_back(std::make_pair(mc[i],       "hist"));
    else{
      histAndOpts.push_back(std::make_pair(mc_witherr[i],       "graphe3"));
      histAndOpts.push_back(std::make_pair(mc[i],       "hist"));
    }
        
    //Do by track breakdown
    if(doTracks){
    /*
      histAndOpts.push_back(std::make_pair(mc_1track[i],       "hist"));
      histAndOpts.push_back(std::make_pair(mc_2track[i],       "hist"));
      histAndOpts.push_back(std::make_pair(data_1track[i],       "hist"));
      histAndOpts.push_back(std::make_pair(data_2track[i],       "hist"));
    */
    }
    
    
    //do by physics process
    else{
      if(doRecoilFit){
	histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_nptune[i], "hist"));
	histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_nntune[i], "hist"));
	histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_qetune[i], "hist"));
	histAndOpts.push_back(std::make_pair(mc_qelike_2p2h[i], "hist"));
      }
      else{
	if(!withresfsi&&!withqefsi){
	  histAndOpts.push_back(std::make_pair(mc_qelike_qe[i],       "hist"));
	  histAndOpts.push_back(std::make_pair(mc_qelike_res[i],       "hist"));
	  histAndOpts.push_back(std::make_pair(mc_qelike_dis[i],       "hist"));
	  histAndOpts.push_back(std::make_pair(mc_qelike_2p2h[i],       "hist"));
	  histAndOpts.push_back(std::make_pair(mc_bkg[i],       "hist"));
	}
	if(withlowrecoilremoved) histAndOpts.push_back(std::make_pair(mc_qelike_2p2h_no_lowrec[i],       "hist"));
	if(withresfsi){
	  histAndOpts.push_back(std::make_pair(mc_qelike_res[i],       "hist"));
	  histAndOpts.push_back(std::make_pair(mc_qelike_res_proton_fsi[i],       "hist"));
	  histAndOpts.push_back(std::make_pair(mc_qelike_res_neutron_fsi[i],       "hist"));
	}
	if(withqefsi){
	  histAndOpts.push_back(std::make_pair(mc_qelike_qe[i],       "hist"));
	  histAndOpts.push_back(std::make_pair(mc_qelike_qe_proton_fsi[i],       "hist"));
	  histAndOpts.push_back(std::make_pair(mc_qelike_qe_neutron_fsi[i],       "hist"));
	}
	

      }
    }
    
    histAndOpts.push_back(std::make_pair(dataStat[i], "histpe"));
    histAndOpts.push_back(std::make_pair(data[i],     "histpe"));
    
    
    
    // ----------------------------------------------------------------------------------
    //
    // First make pt in bins of pz
    
    // Values to multiply each bin by to get them on a similar range

    vector<double> multipliers = GetScales(histAndOpts,false,5,0.75,false);
    GridCanvas* gc=NULL;
    gc = plotYAxis1D(histAndOpts,"P_{t} (GeV/c)","Vis. E (GeV)", doMultipliers ? &multipliers[0] : NULL);
    // Set the y range manually. Can also use gc->Remax() to guess automatically
    if(doRatio) gc->SetYLimits(0,1.99);
    else gc->SetYLimits(0, 5);
    if(doZoom)gc->SetXLimits(0,0.4);
    if(doRatio)gc->SetYTitle("Ratio to MINERvA Tune v1");
    else gc->SetYTitle("Event Rate(x10^{-5}) per GeV^{3}/c^{2}");
    gc->Modified();
    // Example of adding a legend. The co-ordinate system is NDC on the
    // entire canvas, ie (0,0) in the bottom left corner of the canvas
    // (not the individual pad), and (1,1) in the top right
    TLegend* leg=new TLegend(0.7, 0.1, 0.9, 0.3);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(data[i], "MINERvA data", "lpe");
    leg->AddEntry(mc[i], "MINERvA Tune v1", "l");
    if(doTracks){
      /*
      leg->AddEntry(mc_1track[i], "MC 1-track","l");
      leg->AddEntry(mc_2track[i], "MC Multi-track","l");
      leg->AddEntry(data_1track[i], "data 1-track","l");
      leg->AddEntry(data_2track[i], "data Multi-track","l");
      */
    }
    else{
      if(doRecoilFit){
	leg->AddEntry(mc_qelike_2p2h_nptune[i],"2p2h np tune","l");
	leg->AddEntry(mc_qelike_2p2h_nntune[i],"2p2h nn tune","l");
	leg->AddEntry(mc_qelike_2p2h_qetune[i],"2p2h qe tune","l");
	leg->AddEntry(mc_qelike_2p2h[i],"2p2h fraction","l");
      }
      else{
	leg->AddEntry(mc_qelike_qe[i],"QE","l");
	leg->AddEntry(mc_qelike_res[i],"Resonant","l");
	leg->AddEntry(mc_qelike_dis[i],"DIS","l");
	leg->AddEntry(mc_qelike_2p2h[i],"2p2h","l");
	if(withlowrecoilremoved)leg->AddEntry(mc_qelike_2p2h_no_lowrec[i],"2p2h without fit","l");
	if(withresfsi){
	  leg->AddEntry(mc_qelike_res_proton_fsi[i],"RES proton","l");
	  leg->AddEntry(mc_qelike_res_neutron_fsi[i],"RES_neutron","l");	
	}
	if(withqefsi){
	  leg->AddEntry(mc_qelike_qe_proton_fsi[i],"QE proton","l");
	  leg->AddEntry(mc_qelike_qe_neutron_fsi[i],"QE_neutron","l");
	}
	leg->AddEntry(mc_bkg[i],"Background","l");
      }
      
    }

    TLatex mytex;
    mytex.SetTextSize(0.05);
    string mystring =     Form("%.2f < P_{||} (GeV/c) < %.2f",pz3Dbins[i-1],pz3Dbins[i]);
    mytex.DrawLatex(0.35,0.96,mystring.c_str());    
    leg->Draw("SAME");
    if(doTracks){
      if(doBkgSub){
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pt-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pt-bin-%d.eps",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pt-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pt-bin-%d.png",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pt-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pt-bin-%d.C",i));
      }
      else{
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pt-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-trackmult-pt-bin-%d.eps",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pt-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-trackmult-pt-bin-%d.png",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pt-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-trackmult-pt-bin-%d.C",i));
      }
    }
    else{
      if(doZoom){
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d-zoomed.eps",i) : Form("nu-2d-evtrate-model-pt-bin-%d-zoomed.eps",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d-zoomed.png",i) : Form("nu-2d-evtrate-model-pt-bin-%d-zoomed.png",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d-zoomed.C",i) : Form("nu-2d-evtrate-model-pt-bin-%d-zoomed.C",i));
      }
      else{
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-model-pt-bin-%d.eps",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-model-pt-bin-%d.png",i));
	gc->Print(doMultipliers ? Form("nu-2d-evtrate-model-pt-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-model-pt-bin-%d.C",i));
      }
    }
    
    // ----------------------------------------------------------------------------------
    //
    // Now make pz in bins of pt. It's all the same
    
    // Values to multiply each bin by to get them on a similar range
    double  multipliers5[]={4,1,1,1,4,200,75,1,1,1,1,1,1,1};
    
			    
    
    // plotXAxis1D fiddles the x axis values to squash up the tail so it
    // doesn't take up all the horizontal space.
    vector<double> multipliers2 = GetScales(histAndOpts,true,5,0.75,false);
    GridCanvas* gc2=NULL;
    gc2 = plotXAxis1D(histAndOpts,"Visible Energy (GeV)","P_{t} (GeV/c)", doMultipliers ? &multipliers2[0] : NULL);

    mytex.DrawLatex(0.35,0.96,mystring.c_str());
    if(doRatio) gc2->SetYLimits(0,1.99);
    else gc2->SetYLimits(0, 5);
    if(doZoom)gc2->SetXLimits(0,0.4);
    if(doRatio)gc2->SetYTitle("Ratio to MINERvA Tune v1");
    else gc2->SetYTitle("Event Rate(x10^{-5}) per GeV^{3}/c^{2}");

    leg->SetX1(0.77);
    leg->SetY1(0.1);
    leg->SetX2(1.0);
    leg->SetY2(0.4);
    leg->Draw("SAME");
    gc2->Modified();
    if(doTracks){
      if(doBkgSub){
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pz-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pz-bin-%d.eps",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pz-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pz-bin-%d.png",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-bkgsub-trackmult-pz-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-bkgsub-trackmult-pz-bin-%d.C",i));
      }
      else{
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pz-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-trackmult-pz-bin-%d.eps",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pz-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-trackmult-pz-bin-%d.png",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-trackmult-pz-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-trackmult-pz-bin-%d.C",i));
      }
    }
    else{
      if(doZoom){
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d-zoomed.eps",i) : Form("nu-2d-evtrate-model-pz-bin-%d-zoomed.eps",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d-zoomed.png",i) : Form("nu-2d-evtrate-model-pz-bin-%d-zoomed.png",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d-zoomed.C",i) : Form("nu-2d-evtrate-model-pz-bin-%d-zoomed.C",i));
      }
      else{
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.eps",i) : Form("nu-2d-evtrate-model-pz-bin-%d.eps",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.png",i) : Form("nu-2d-evtrate-model-pz-bin-%d.png",i));
	gc2->Print(doMultipliers ? Form("nu-2d-evtrate-model-pz-multiplier-bin-%d.C",i) : Form("nu-2d-evtrate-model-pz-bin-%d.C",i));
      }
    }
  }
}

int main(int argc, char* argv[])
{

  string location = argv[1];
  string s_withlowrecoilremoved = argv[2];
  string s_withresfsi = argv[3];
  string s_withqefsi = argv[4];
  string s_dorecoilfit = argv[5];
  string s_ratio = argv[6];
  string s_zoom = argv[7];
  string s_mcerr = argv[8];
  string s_bkgconst = argv[9];

  bool withlowrecoilremoved = s_withlowrecoilremoved=="1" ? true:false;
  bool withresfsi = s_withresfsi=="1" ? true:false;
  bool withqefsi = s_withqefsi=="1" ? true:false;
  bool doRecoilFit = s_dorecoilfit=="1" ? true:false;
  bool doRatio = s_ratio=="1" ? true:false;
  bool doZoom = s_zoom=="1"? true:false;
  bool domcerr = s_mcerr=="1"? true:false;
  bool doconst = s_bkgconst=="1"?true:false;
  //multipliers
  //  makePlots(true,true,false,location,withlowrecoilremoved,withresfsi,withqefsi,doRecoilFit,false);
  //  makePlots(true,true,true,location,withlowrecoilremoved,withresfsi,withqefsi,doRecoilFit,false);
  makePlots(true,false,false,location,withlowrecoilremoved,withresfsi,withqefsi,doRecoilFit,false,doZoom,domcerr,doconst);
  //standard
  //  makePlots(false,true,false,location,withlowrecoilremoved,withresfsi,withqefsi,doRecoilFit,doRatio);
  //  makePlots(false,true,true,location,withlowrecoilremoved,withresfsi,withqefsi,doRecoilFit,doRatio);
  makePlots(false,false,false,location,withlowrecoilremoved,withresfsi,withqefsi,doRecoilFit,doRatio,doZoom,domcerr,doconst);

  return 0;
}
