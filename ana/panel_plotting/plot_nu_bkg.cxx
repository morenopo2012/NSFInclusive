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
#include <iostream>
#include <fstream>
#include <string>

using namespace PlotUtils;
void Create3DBkgConst(MnvH2D*&h, MnvH2D* h_temp, HyperDimLinearizer *hd, vector<double> zcoords){
  
  //Get bin centers z
  vector<double> zcen;
  for(uint i=0;i<zcoords.size();i++){
    if(i==0)zcen.push_back(1.0); //underflow
    if(i==zcoords.size()-1) zcen.push_back(zcen.back()+20); //overflow
    else zcen.push_back((zcoords[i+1]+zcoords[i])/2); //normal
    


  }
  //  for(uint i=0;i<zcen.size();i++)    cout << i << "\t" << zcen[i]<<endl;
  //vert and lat errors
  vector<string> verrors = h_temp->GetVertErrorBandNames();
  vector<string> lerrors = h_temp->GetLatErrorBandNames();

  //get pt and recoil coords
  int num_recoil_bins = h_temp->GetNbinsX()+2;
  int num_pt_bins = h_temp->GetNbinsY()+2;
  
  //CV
  for(uint z=0;z<zcen.size();z++){//loop over Z
    for(int i=0;i<num_recoil_bins;i++){//recoil
      for(int j=0;j<num_pt_bins;j++){//pt
	//Get 3D coord
	double scale = h_temp->GetBinContent(i,j);
	if(scale==0) scale=1;
	vector<double> vals;
	vals.push_back(h_temp->GetXaxis()->GetBinCenter(i)/1000.0);
	vals.push_back(h_temp->GetYaxis()->GetBinCenter(j));
	vals.push_back(zcen[z]);
	
	pair<int,int> target = hd->GetBin(vals);
	//	cout << i << "\t" << j << "\t" << z << "\t" <<h_temp->GetXaxis()->GetBinCenter(i) << "\t" << h_temp->GetYaxis()->GetBinCenter(j) << "\t" << zcen[z] << "\t" << target.first << "\t" << target.second <<"\t" <<h_temp->GetBinContent(i,j)<< endl;
	h->SetBinContent(target.first+1,target.second,scale);
	h->SetBinError(target.first,j,h_temp->GetBinError(i,j));

      }
    }
  }

  //Populate errorsbands and the CV of each
  h->AddMissingErrorBandsAndFillWithCV(*h_temp);

  //loop over verts
  for(uint v=0;v<verrors.size();v++){
    MnvVertErrorBand2D *vband_temp = h_temp->GetVertErrorBand(verrors[v]);
    MnvVertErrorBand2D *vband = h->GetVertErrorBand(verrors[v]);
    int nuni = vband_temp->GetNHists();
    for(int u=0;u<nuni;u++){
      TH2D *h_vband_temp = vband_temp->GetHist(u);
      TH2D *h_vband = vband->GetHist(u);
      
      for(uint z=0;z<zcen.size();z++){//loop over Z
	for(int i=0;i<num_recoil_bins;i++){//recoil
	  for(int j=0;j<num_pt_bins;j++){//pt
	    //Get 3D coord
	    double scale = h_vband_temp->GetBinContent(i,j);
	    if(scale==0) scale=1;
	    vector<double> vals;
	    vals.push_back(h_vband_temp->GetXaxis()->GetBinCenter(i)/1000.0);
	    vals.push_back(h_vband_temp->GetYaxis()->GetBinCenter(j));
	    vals.push_back(zcen[z]);
	    
	    pair<int,int> target = hd->GetBin(vals);
	    h_vband->SetBinContent(target.first+1,target.second,scale);
	  }
	}
      }
    }
  }

  //loop over lats
  for(uint l=0;l<lerrors.size();l++){
    MnvLatErrorBand2D *lband_temp = h_temp->GetLatErrorBand(lerrors[l]);
    MnvLatErrorBand2D *lband = h->GetLatErrorBand(lerrors[l]);
    int nuni = lband_temp->GetNHists();
    for(int u=0;u<nuni;u++){
      TH2D *h_lband_temp = lband_temp->GetHist(u);
      TH2D *h_lband = lband->GetHist(u);
      
      for(uint z=0;z<zcen.size();z++){//loop over Z
	for(int i=0;i<num_recoil_bins;i++){//recoil
	  for(int j=0;j<num_pt_bins;j++){//pt
	    //Get 3D coord
	    vector<double> vals;
	    double scale = h_lband_temp->GetBinContent(i,j);
	    if(scale==0) scale=1;
	    vals.push_back(h_lband_temp->GetXaxis()->GetBinCenter(i)/1000.0);
	    vals.push_back(h_lband_temp->GetYaxis()->GetBinCenter(j));
	    vals.push_back(zcen[z]);
	    
	    pair<int,int> target = hd->GetBin(vals);
	    h_lband->SetBinContent(target.first+1,target.second,scale);
	  }
	}
      }
    }
  }
}


void makePlots(bool doMultipliers, string location,string sample, bool doRatio, bool constrainSig)
{
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
  pz3Dbins.push_back(8.0);
  //These are added to fix unfolding
  pz3Dbins.push_back(10.0);//fix unf
  pz3Dbins.push_back(20.0);

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



  TFile f1(Form("%s/SideBandFit_SVD.root",location.c_str()));//scale
  TFile f2(Form("%s/MuonEventSelection_MakeFlux-1_Multiplicity-0_Sample-%s_CombinedPlaylists.root",location.c_str(),sample.c_str()));//sample 1 trk
  //TFile f3(Form("%s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-%s_CombinedPlaylists.root",location.c_str(),sample.c_str()));//sample 2 trk


  string scaleform = "Signal";
  if(sample=="MichelSideBand") scaleform = "Michel";
  if(sample=="BlobSideBand") scaleform = "Blob";
  if(sample=="MicBlobSideBand") scaleform = "Micblob";
  
  MnvH2D* scale_1track_temp = (MnvH2D*)f1.Get(Form("%s_qelikenot_scaled",scaleform.c_str()));
  MnvH2D* scale_2track_temp = (MnvH2D*)f1.Get(Form("%s_qelikenot_scaled",scaleform.c_str()));

  MnvH2D* scale_1track_sigtemp = (MnvH2D*)f1.Get("qescale");
  MnvH2D* scale_2track_sigtemp = (MnvH2D*)f1.Get("qescale");



  MnvH2D* dataMnv_1trk=(MnvH2D*)f2.Get("h_pzptrec_data");
  //  MnvH2D* dataMnv_2trk=(MnvH2D*)f3.Get("h_pzptrec_data");
  
  MnvH2D* mcsigMnv_1trk=(MnvH2D*)f2.Get("h_pzptrec_qelike");
  //  MnvH2D* mcsigMnv_2trk=(MnvH2D*)f3.Get("h_pzptrec_qelike");

  MnvH2D* mcbkgMnv_1trk=(MnvH2D*)f2.Get("h_pzptrec_qelikenot");
  //  MnvH2D* mcbkgMnv_2trk=(MnvH2D*)f3.Get("h_pzptrec_qelikenot");
  /*
  MnvH2D* mcbkgMnv_scp_1trk=(MnvH2D*)f2.Get("h_pzptrec_qelikenot_singlechargedpion");
  MnvH2D* mcbkgMnv_scp_2trk=(MnvH2D*)f3.Get("h_pzptrec_qelikenot_singlechargedpion");

  MnvH2D* mcbkgMnv_snp_1trk=(MnvH2D*)f2.Get("h_pzptrec_qelikenot_singleneutralpion");
  MnvH2D* mcbkgMnv_snp_2trk=(MnvH2D*)f3.Get("h_pzptrec_qelikenot_singleneutralpion");

  MnvH2D* mcbkgMnv_mp_1trk=(MnvH2D*)f2.Get("h_pzptrec_qelikenot_multipion");
  MnvH2D* mcbkgMnv_mp_2trk=(MnvH2D*)f3.Get("h_pzptrec_qelikenot_multipion");

  MnvH2D* mcbkgMnv_nop_1trk=(MnvH2D*)f2.Get("h_pzptrec_qelikenot_not_scp_snp_mp");
  MnvH2D* mcbkgMnv_nop_2trk=(MnvH2D*)f3.Get("h_pzptrec_qelikenot_not_scp_snp_mp");
  */
  //unscaled
  MnvH2D* mcMnv_1trk=(MnvH2D*)f2.Get("h_pzptrec_mc");
  //  MnvH2D* mcMnv_2trk=(MnvH2D*)f3.Get("h_pzptrec_mc");

  MnvH2D* scale_1track = dataMnv_1trk->Clone("1trackscale");
  scale_1track->ClearAllErrorBands();
  scale_1track->Reset();
  //  MnvH2D*scale_2track = dataMnv_2trk->Clone("2trackscale");

  MnvH2D* scale_sig_1track = dataMnv_1trk->Clone("sig1trackscale");
  //  MnvH2D*scale_sig_2track = dataMnv_2trk->Clone("sig2trackscale");

  cout << "Make 3D constraint" << endl;
  Create3DBkgConst(scale_1track, scale_1track_temp, my3d, pz3Dbins );
  cout << "done1"<< endl;
  // Create3DBkgConst(scale_2track, scale_2track_temp, my3d, pz3Dbins );
  cout << "Done" << endl;

  cout << "Make 3D constraint" << endl;
  Create3DBkgConst(scale_sig_1track, scale_1track_sigtemp, my3d, pz3Dbins );
  cout << "done1"<< endl;
  //  Create3DBkgConst(scale_sig_2track, scale_2track_sigtemp, my3d, pz3Dbins );
  cout << "Done" << endl;

  //apply scale
  mcbkgMnv_1trk->Multiply(scale_1track,mcbkgMnv_1trk);
  //  mcbkgMnv_2trk->Multiply(scale_2track,mcbkgMnv_2trk);
  /*
  mcbkgMnv_scp_1trk->Multiply(scale_1track,mcbkgMnv_scp_1trk);
  mcbkgMnv_scp_2trk->Multiply(scale_2track,mcbkgMnv_scp_2trk);
  mcbkgMnv_snp_1trk->Multiply(scale_1track,mcbkgMnv_snp_1trk);
  mcbkgMnv_snp_2trk->Multiply(scale_2track,mcbkgMnv_snp_2trk);
  mcbkgMnv_mp_1trk->Multiply(scale_1track,mcbkgMnv_mp_1trk);
  mcbkgMnv_mp_2trk->Multiply(scale_2track,mcbkgMnv_mp_2trk);
  mcbkgMnv_nop_1trk->Multiply(scale_1track,mcbkgMnv_nop_1trk);
  mcbkgMnv_nop_2trk->Multiply(scale_2track,mcbkgMnv_nop_2trk);
  */
  //add samples (total rate)
  if(!constrainSig){//All data All Bkg
    mcsigMnv_1trk->Add(mcbkgMnv_1trk);
    //    mcsigMnv_2trk->Add(mcbkgMnv_2trk);
  }
  else{//MC - Sig - "other"
    mcsigMnv_1trk->Multiply(mcsigMnv_1trk,scale_sig_1track);
    //    mcsigMnv_2trk->Multiply(mcsigMnv_2trk,scale_sig_2track);

    mcsigMnv_1trk->Add(mcbkgMnv_1trk);
    //    mcsigMnv_2trk->Add(mcbkgMnv_2trk);

  }
  
  //Background fraction
  //  mcbkgMnv_1trk->Divide(mcbkgMnv_1trk,mcsigMnv_1trk);
  //  mcbkgMnv_2trk->Divide(mcbkgMnv_2trk,mcsigMnv_2trk);
  




  MnvH2D *dataMnv_alltrk = dataMnv_1trk->Clone("data_all");
  MnvH2D *mcMnv_alltrk = mcMnv_1trk->Clone("all_unconstrained");
  MnvH2D *mcsigMnv_alltrk = mcsigMnv_1trk->Clone("all_constrained");

  //  dataMnv_alltrk->Add(dataMnv_2trk);
  //  mcMnv_alltrk->Add(mcMnv_2trk);
  //  mcsigMnv_alltrk->Add(mcsigMnv_2trk);

  TFile *out = new TFile(Form("ScaledBkgCat_%s.root",sample.c_str()),"RECREATE");
  mcsigMnv_1trk->Write("h_pzptrec_mc_1trk");
  //  mcsigMnv_2trk->Write("h_pzptrec_mc_2trk");
  mcsigMnv_alltrk->Write("h_pzptrec_mc_alltrk");
  dataMnv_1trk->Write("h_pzptrec_data_1trk");
  //  dataMnv_2trk->Write("h_pzptrec_data_2trk");
  dataMnv_alltrk->Write("h_pzptrec_data_alltrk");
  mcbkgMnv_1trk->Write("h_pzptrec_mc_1trk_bkgfraction");
  //  mcbkgMnv_2trk->Write("h_pzptrec_mc_2trk_bkgfraction");
  mcMnv_alltrk->Write("h_pzptrec_mc_alltrk_unconstrained");
  scale_1track->Write("Scale");
  scale_1track_temp->Write("OrigScale");
  out->Close();


  std::cout << "Starting up getting the projections" << std::endl;
  /*
  std::vector<TH2D*> dataresults_1trk = my3d->Get2DHistos(dataMnv_1trk,true);
  std::vector<TH2D*> mcresults_1trk = my3d->Get2DHistos(mcMnv_1trk,true);
  std::vector<TH2D*> scaledmcresults_1trk = my3d->Get2DHistos(mcsigMnv_1trk,true);

  std::vector<TH2D*> dataresults_2trk = my3d->Get2DHistos(dataMnv_2trk,true);
  std::vector<TH2D*> mcresults_2trk = my3d->Get2DHistos(mcMnv_2trk,true);
  std::vector<TH2D*> scaledmcresults_2trk = my3d->Get2DHistos(mcsigMnv_2trk,true);
  */
  std::vector<TH2D*> dataresults_alltrk = my3d->Get2DHistos(dataMnv_alltrk,true);
  std::vector<TH2D*> mcresults_alltrk = my3d->Get2DHistos(mcMnv_alltrk,true);
  std::vector<TH2D*> scaledmcresults_alltrk = my3d->Get2DHistos(mcsigMnv_alltrk,true);


  
  for(unsigned int i=1;i<pz3Dbins.size();i++){
    double pzwidth = pz3Dbins[i]-pz3Dbins[i-1];
    /*
    dataresults_1trk[i]->Scale(1/pzwidth,"width");
    mcresults_1trk[i]->Scale(1/pzwidth,"width");
    scaledmcresults_1trk[i]->Scale(1/pzwidth,"width");

    dataresults_2trk[i]->Scale(1/pzwidth,"width");
    mcresults_2trk[i]->Scale(1/pzwidth,"width");
    scaledmcresults_2trk[i]->Scale(1/pzwidth,"width");
    */
    dataresults_alltrk[i]->Scale(1/pzwidth,"width");
    mcresults_alltrk[i]->Scale(1/pzwidth,"width");
    scaledmcresults_alltrk[i]->Scale(1/pzwidth,"width");
  }
  
  vector<int> mycolors = getColors(2);
  for(unsigned int i=1;i<pz3Dbins.size();i++){
    cout << "DOING BIN " << i << endl;
    /*
    mcresults_1trk[i]->SetLineColor(kRed);
    mcresults_1trk[i]->SetLineWidth(2);
    mcresults_2trk[i]->SetLineColor(kRed);
    mcresults_2trk[i]->SetLineWidth(2);
    */
    mcresults_alltrk[i]->SetLineColor(kRed);
    mcresults_alltrk[i]->SetLineWidth(2);

    //mcresults_1trk[i]->SetFillStyle(3005);
    //mcresults_2trk[i]->SetFillStyle(3005);
    mcresults_alltrk[i]->SetFillStyle(3005);

    //mcresults_1trk[i]->SetFillColor(kRed);
    //mcresults_2trk[i]->SetFillColor(kRed);
    mcresults_alltrk[i]->SetFillColor(kRed);
    /*
    scaledmcresults_1trk[i]->SetLineColor(kBlue);
    scaledmcresults_1trk[i]->SetLineWidth(2);
    scaledmcresults_2trk[i]->SetLineColor(kBlue);
    scaledmcresults_2trk[i]->SetLineWidth(2);
    */
    scaledmcresults_alltrk[i]->SetLineColor(kBlue);
    scaledmcresults_alltrk[i]->SetLineWidth(2);

    //scaledmcresults_1trk[i]->SetFillStyle(3004);
    //scaledmcresults_2trk[i]->SetFillStyle(3004);
    scaledmcresults_alltrk[i]->SetFillStyle(3004);


    //scaledmcresults_1trk[i]->SetFillColor(kBlue);
    //scaledmcresults_2trk[i]->SetFillColor(kBlue);
    scaledmcresults_alltrk[i]->SetFillColor(kBlue);

    /*
    dataresults_1trk[i]->SetMarkerStyle(kFullCircle);
    dataresults_1trk[i]->SetMarkerSize(0.7);
    dataresults_1trk[i]->SetLineColor(kBlack);
    dataresults_1trk[i]->SetLineWidth(2);
    dataresults_2trk[i]->SetMarkerStyle(kFullCircle);
    dataresults_2trk[i]->SetMarkerSize(0.7);
    dataresults_2trk[i]->SetLineColor(kBlack);
    dataresults_2trk[i]->SetLineWidth(2);
    */
    dataresults_alltrk[i]->SetMarkerStyle(kFullCircle);
    dataresults_alltrk[i]->SetMarkerSize(0.7);
    dataresults_alltrk[i]->SetLineColor(kBlack);
    dataresults_alltrk[i]->SetLineWidth(2);

    /*
    if(true){
      for(int j=1;j<7;j++){
	double dr1 = dataresults_1trk[i]->ProjectionX("tmp1",j,j)->Integral();
	double mr1 = mcresults_1trk[i]->ProjectionX("tmp2",j,j)->Integral();
	double smr1 = scaledmcresults_1trk[i]->ProjectionX("tmp1",j,j)->Integral();

	double dr2 = dataresults_2trk[i]->ProjectionX("tmp1",j,j)->Integral();
	double mr2 = mcresults_2trk[i]->ProjectionX("tmp2",j,j)->Integral();
	double smr2 = scaledmcresults_2trk[i]->ProjectionX("tmp1",j,j)->Integral();

	cout << "****" << endl;
	cout << j << "\t" << dr1  << "\t" << mr1<< "\t" << smr1 << "\t" << smr1/mr1 << "\t" <<scale_1track->GetBinContent(1,j) << endl;
	cout << j << "\t" << dr2  << "\t" << mr2<< "\t" << smr2 << "\t" << smr2/mr2 << "\t" <<scale_2track->GetBinContent(1,j) << endl;
	cout << "****" << endl;
      }
    }
    */


    if(doRatio){//scale eventhing to scaled bkg

      //      TH2D *tmp_1trk = (TH2D*)scaledmcresults_1trk[i]->Clone();
      //      TH2D *tmp_2trk = (TH2D*)scaledmcresults_2trk[i]->Clone();
      TH2D *tmp_alltrk = (TH2D*)scaledmcresults_alltrk[i]->Clone();
      
      int xbins = tmp_alltrk->GetNbinsX()+2;
      int ybins = tmp_alltrk->GetNbinsY()+2;
      
      for(int i=0;i<xbins;i++){
	for(int j=0;j<ybins;j++){
	  //	  tmp_1trk->SetBinError(i,j,0);
	  //	  tmp_2trk->SetBinError(i,j,0);
	  tmp_alltrk->SetBinError(i,j,0);
	}
      }
      
      /*
      dataresults_1trk[i]->Divide(tmp_1trk);
      mcresults_1trk[i]->Divide(tmp_1trk);
      scaledmcresults_1trk[i]->Divide(tmp_1trk);
      
      dataresults_2trk[i]->Divide(tmp_2trk);
      mcresults_2trk[i]->Divide(tmp_2trk);
      scaledmcresults_2trk[i]->Divide(tmp_2trk);
      */
      dataresults_alltrk[i]->Divide(tmp_alltrk);
      mcresults_alltrk[i]->Divide(tmp_alltrk);
      scaledmcresults_alltrk[i]->Divide(tmp_alltrk);
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
    /*
    std::vector<std::pair<TH2*, const char*> > histAndOpts_1trk;
    histAndOpts_1trk.push_back(std::make_pair(dataresults_1trk[i], "histpe"));
    histAndOpts_1trk.push_back(std::make_pair(mcresults_1trk[i],       "histle3"));
    histAndOpts_1trk.push_back(std::make_pair(scaledmcresults_1trk[i],       "histle3"));


    std::vector<std::pair<TH2*, const char*> > histAndOpts_2trk;
    histAndOpts_2trk.push_back(std::make_pair(dataresults_2trk[i], "histpe"));
    histAndOpts_2trk.push_back(std::make_pair(mcresults_2trk[i],       "histle3"));
    histAndOpts_2trk.push_back(std::make_pair(scaledmcresults_2trk[i],       "histle3"));

    */
    std::vector<std::pair<TH2*, const char*> > histAndOpts_alltrk;
    histAndOpts_alltrk.push_back(std::make_pair(dataresults_alltrk[i], "histpe"));
    histAndOpts_alltrk.push_back(std::make_pair(mcresults_alltrk[i],       "histle3"));
    histAndOpts_alltrk.push_back(std::make_pair(scaledmcresults_alltrk[i],       "histle3"));


    cout << "Mults" << endl;


    // ----------------------------------------------------------------------------------
    //
    // First make pt in bins of pz
      
    



    //    vector<double> multi_1 = GetScales(histAndOpts_1trk, true, 100000,0.75);

    // Example of adding a legend. The co-ordinate system is NDC on the
    // entire canvas, ie (0,0) in the bottom left corner of the canvas
    // (not the individual pad), and (1,1) in the top right
    TLegend* leg=new TLegend(0.79, 0.1, 1, 0.4);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(dataresults_alltrk[i], "MINERvA data", "lpe");
    leg->AddEntry(mcresults_alltrk[i], "Minerva Tune v1", "l");
    leg->AddEntry(scaledmcresults_alltrk[i], "Tune v1 w/ sc. bkg.", "l");
    

    TLatex mytex;
    mytex.SetTextSize(0.05);
    string mystring =     Form("%.2f < P_{||} (GeV/c) < %.2f",pz3Dbins[i-1],pz3Dbins[i]);
    /*
    GridCanvas* gc1=NULL;
    gc1=plotXAxis1D(histAndOpts_1trk, "Visible Energy (GeV)","P_{t}", doMultipliers ? &multi_1 [0]: NULL);
    mytex.DrawLatex(0.35,0.96,mystring.c_str());

    // Set the y range manually. Can also use gc1->Remax() to guess automatically
    if(doRatio&&!constrainSig)gc1->SetYLimits(0,1.99);
    else if(doRatio&&constrainSig) gc1->SetYLimits(0,1.99);
    else gc1->SetYLimits(0, 100000);
    //gc1->SetYLimits(0, 1e-39);
    if(doRatio)gc1->SetYTitle("Ratio to Minerva Tune v1");
    else gc1->SetYTitle("Events per GeV^{3}/c^{2}");
    leg->Draw("SAME");
    gc1->Modified();
  

    if(doRatio){
      if(constrainSig){
	gc1->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-sigconstraint-ratio-1track.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-ratio-1track.eps",sample.c_str(),i));
	gc1->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-ratio-1track.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-ratio-1track.png",sample.c_str(),i));
	gc1->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-ratio-1track.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-ratio-1track.C",sample.c_str(),i));
      }
      else{
	gc1->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-ratio-1track.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-ratio-1track.eps",sample.c_str(),i));
	gc1->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-ratio-1track.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-ratio-1track.png",sample.c_str(),i));
	gc1->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-ratio-1track.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-ratio-1track.C",sample.c_str(),i));
      }
    }
    else{
      if(constrainSig){
	gc1->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-sigconstraint-1track.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-1track.eps",sample.c_str(),i));
	gc1->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-1track.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-1track.png",sample.c_str(),i));
	gc1->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-1track.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-1track.C",sample.c_str(),i));
      }
      else{
	gc1->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-1track.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-1track.eps",sample.c_str(),i));
	gc1->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-1track.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-1track.png",sample.c_str(),i));
	gc1->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-1track.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-1track.C",sample.c_str(),i));
      }
    }

    



    vector<double> multi_2 = GetScales(histAndOpts_2trk, true,100000,0.75);
    GridCanvas* gc2=NULL;
    gc2=plotXAxis1D(histAndOpts_2trk, "Visible Energy (GeV)","P_{t}", doMultipliers ? &multi_2[0] : NULL);
    mytex.DrawLatex(0.35,0.96,mystring.c_str());

    // Set the y range manually. Can also use gc2->Remax() to guess automatically
    if(doRatio&&!constrainSig)gc2->SetYLimits(0,1.99);
    else if(doRatio&&constrainSig) gc2->SetYLimits(0,1.99);
    else gc2->SetYLimits(0, 100000);
    //gc2->SetYLimits(0, 1e-39);

    if(doRatio)gc2->SetYTitle("Ratio to Minerva Tune v1");
    else gc2->SetYTitle("Events per GeV^{3}/c^{2}");
    leg->Draw("SAME");
    gc2->Modified();
  
    if(doRatio){
      if(constrainSig){
	gc2->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-sigconstraint-ratio-2track.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-ratio-2track.eps",sample.c_str(),i));
	gc2->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-ratio-2track.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-ratio-2track.png",sample.c_str(),i));
	gc2->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-ratio-2track.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-ratio-2track.C",sample.c_str(),i));
      }
      else{
	gc2->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-ratio-2track.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-ratio-2track.eps",sample.c_str(),i));
	gc2->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-ratio-2track.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-ratio-2track.png",sample.c_str(),i));
	gc2->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-ratio-2track.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-ratio-2track.C",sample.c_str(),i));
      }
    }
    else{
      if(constrainSig){
	gc2->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-sigconstraint-2track.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-2track.eps",sample.c_str(),i));
	gc2->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-2track.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-2track.png",sample.c_str(),i));
	gc2->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-2track.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-2track.C",sample.c_str(),i));
      }
      else{
	gc2->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-2track.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-2track.eps",sample.c_str(),i));
	gc2->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-2track.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-2track.png",sample.c_str(),i));
	gc2->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-2track.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-2track.C",sample.c_str(),i));
      }
    }
    */

    vector<double> multi_3 = GetScales(histAndOpts_alltrk, true, 100000,0.75);
    GridCanvas* gc3=NULL;
    gc3=plotXAxis1D(histAndOpts_alltrk, "Visible Energy (GeV)","P_{t}", doMultipliers ?  &multi_3[0] : NULL);


    mytex.DrawLatex(0.35,0.96,mystring.c_str());

    // Set the y range manually. Can also use gc3->Remax() to guess automatically
    if(doRatio&&!constrainSig)gc3->SetYLimits(0,1.99);
    else if(doRatio&&constrainSig)gc3->SetYLimits(0,1.99);
    else gc3->SetYLimits(0, 100000);
    //gc3->SetYLimits(0, 1e-39);

    if(doRatio)gc3->SetYTitle("Ratio to scaled Minerva Tune v1");
    else gc3->SetYTitle("Events per MeV#timesGeV^{2}/c^{2}");
    leg->Draw("SAME");
    gc3->Modified();
  
    if(doRatio){
      if(constrainSig){
	gc3->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-sigconstraint-ratio-alltrack.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-ratio-alltrack.eps",sample.c_str(),i));
	gc3->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-ratio-alltrack.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-ratio-alltrack.png",sample.c_str(),i));
	gc3->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-ratio-alltrack.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-ratio-alltrack.C",sample.c_str(),i));
      }
      else{
	gc3->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-ratio-alltrack.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-ratio-alltrack.eps",sample.c_str(),i));
	gc3->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-ratio-alltrack.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-ratio-alltrack.png",sample.c_str(),i));
	gc3->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-ratio-alltrack.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-ratio-alltrack.C",sample.c_str(),i));
      }
    }
    else{
      if(constrainSig){
	gc3->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-sigconstraint-alltrack.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-alltrack.eps",sample.c_str(),i));
	gc3->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-alltrack.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-alltrack.png",sample.c_str(),i));
	gc3->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-sigconstraint-alltrack.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-sigconstraint-alltrack.C",sample.c_str(),i));
      }
      else{
	gc3->Print(doMultipliers ? Form("nu-3d-bkg_%s-comps-pt-multiplier_bin_%d-alltrack.eps",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-alltrack.eps",sample.c_str(),i));
	gc3->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-alltrack.png",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-alltrack.png",sample.c_str(),i));
	gc3->Print(doMultipliers ? Form("nu-3d-bkg-%s-comps-pt-multiplier_bin_%d-alltrack.C",sample.c_str(),i) : Form("nu-3d-bkg-%s-comps-pt_bin_%d-alltrack.C",sample.c_str(),i));
      }
    }
  }
}

int main(int argc, char* argv[])
{
  
  string s_doRatio = argv[2];
  string s_doSig = argv[3];
  bool doRatio = s_doRatio=="1" ? true:false;
  bool doSigConst = s_doSig=="1" ? true:false;
  if(doRatio){
    makePlots(false,argv[1],"MichelSideBand",doRatio,doSigConst);
    makePlots(false,argv[1],"BlobSideBand",doRatio,doSigConst);
    makePlots(false,argv[1],"MicBlobSideBand",doRatio,doSigConst);
    makePlots(false,argv[1],"Signal",doRatio,doSigConst);
  }
  else{
    makePlots(true,argv[1],"MichelSideBand",doRatio,doSigConst);
    makePlots(true,argv[1],"BlobSideBand",doRatio,doSigConst);
    makePlots(true,argv[1],"MicBlobSideBand",doRatio,doSigConst);
    makePlots(true,argv[1],"Signal",doRatio,doSigConst);
  }

  return 0;
}
