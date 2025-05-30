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
#include "TH3D.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;

void fillExternal(MnvH2D *target, TH3D* source, HyperDimLinearizer* mapper, int bpt, int bpz, int brecoil){
  target->Reset();

  // x=recoil, y=pt, z= pz
  for(int i=0;i<bpt-1;i++){//pt
    for(int j=0;j<bpz-1;j++){//pz
      for(int k=0;k<brecoil-1;k++){//recoil
	double pt = source->GetYaxis()->GetBinCenter(i+1);
	double pz = source->GetZaxis()->GetBinCenter(j+1);
	double recoil = source->GetXaxis()->GetBinCenter(k+1);//Convert MeV to GeV
	if(recoil >1.0) recoil = 1.5;
	double xsec = source->GetBinContent(k+1,i+1,j+1);
	double xsec_err = source->GetBinError(k+1,i+1,j+1);
	

	vector<double> val_vect;
	val_vect.push_back(recoil);
	val_vect.push_back(pt);
	val_vect.push_back(pz);
	float globx_long = mapper->GetBin(val_vect).first+0.0001;
	float globy_long = mapper->GetBin(val_vect).second;
	cout << i << "\t" << j << "\t" << k << "\t" << pt << "\t" << pz << "\t" << recoil << "\t" << globx_long << "\t" << globy_long << "\t" << xsec << "\t" << xsec_err<< endl;
	target->Fill(globx_long,pt,xsec);
	target->SetBinError(globx_long,globy_long,xsec_err);
	
      }
    }
  }
}









void makePlots(bool doMultipliers, string location,bool doLog, bool doRatio, int set, bool doZoom)
{

  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  
  TFile f1(Form("%s_CV_Zexp_Bubble/CrossSection_per_nucleon_3D_pzptreco_iterations_10_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  //  TFile f2(Form("%s_CV_DIS_AMU/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  //  TFile f3(Form("%s_CV_DIS_NCTEQ/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  //  TFile f4(Form("%s_CV_Elastic_FSI/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  //  TFile f5(Form("%s_CV_MK/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  TFile f6(Form("%s_CV_RPA_Res_MINOS/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  //  TFile f7(Form("%s_CV_RPA_Res_Nieves/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  //  TFile f8(Form("%s_CV_ZExpansion/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  TFile f9(Form("%s_CV_minerva_joint_lowq2/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  //  TFile f10(Form("%s_default/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  //  TFile f11(Form("%s_pion_2p2h/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  TFile f12(Form("%s_pion_rpa/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final
  //  TFile f13(Form("%s_piontune/CrossSection_per_nucleon_3D_pzptreco_iterations_1_CombinedPlaylists.root_big3d.root",location.c_str()));//final

  
  //External
  //NEUT
  TFile f14("/pnfs/minerva/persistent/users/drut1186/3D_ME_Results_PostWineCheese/PubAttempt_v26_6Pz_11Pt_13Recoil_NorminalFSIUnc_NeutronCV/CCQELike_Models/Neut_v5_SF_ma103.root");//NEUT 5.4.0 SF is MAQE=1.03
  //  TFile f15("/minerva/data/users/drut1186/3D_MODEL_PREDICTIONS/Oct212019/3d/MINERvA_ME_2_NEUT533_Flat.root_mcstats_mode_output.root");//NEUT 5.3.3 are RFG LS MAQE=1.20
  //  TFile f16("/minerva/data/users/drut1186/3D_MODEL_PREDICTIONS/Oct212019/3d/MINERvA_ME_2_NEUT533_Flat_wRPA.root_mcstats_mode_output.root");//NEUT 5.3.3 are RFG LS MAQE=1.20
  TFile f17("/pnfs/minerva/persistent/users/drut1186/3D_ME_Results_PostWineCheese/PubAttempt_v26_6Pz_11Pt_13Recoil_NorminalFSIUnc_NeutronCV/CCQELike_Models/Neut_v5_LFG_ma105.root");//NEUT 5.4.0 LFG Nieves 1p1h, Nieves 2p2h is MAQE=1.05

  //GENIE
  TFile f18("/pnfs/minerva/persistent/users/drut1186/3D_ME_Results_PostWineCheese/PubAttempt_v26_6Pz_11Pt_13Recoil_NorminalFSIUnc_NeutronCV/CCQELike_Models/GENIE_v3_02a.root");
  TFile f19("/pnfs/minerva/persistent/users/drut1186/3D_ME_Results_PostWineCheese/PubAttempt_v26_6Pz_11Pt_13Recoil_NorminalFSIUnc_NeutronCV/CCQELike_Models/GENIE_v3_02b.root");
  TFile f20("/pnfs/minerva/persistent/users/drut1186/3D_ME_Results_PostWineCheese/PubAttempt_v26_6Pz_11Pt_13Recoil_NorminalFSIUnc_NeutronCV/CCQELike_Models/GENIE_v3_10a.root");
  TFile f21("/pnfs/minerva/persistent/users/drut1186/3D_ME_Results_PostWineCheese/PubAttempt_v26_6Pz_11Pt_13Recoil_NorminalFSIUnc_NeutronCV/CCQELike_Models/GENIE_v3_10b.root");
  //NuWro
  TFile f22("/pnfs/minerva/persistent/users/drut1186/3D_ME_Results_PostWineCheese/PubAttempt_v26_6Pz_11Pt_13Recoil_NorminalFSIUnc_NeutronCV/CCQELike_Models/NuWro_SF_v1902.root");
  TFile f23("/pnfs/minerva/persistent/users/drut1186/3D_ME_Results_PostWineCheese/PubAttempt_v26_6Pz_11Pt_13Recoil_NorminalFSIUnc_NeutronCV/CCQELike_Models/NuWro_LFG_v1902.root");
  //GiBUU
  TFile f24("");

  

  MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_pzptrec_data_nobck_unfold_effcor_cross_section");
  MnvH2D* mcMnv=(MnvH2D*)f1.Get("h_pzptrec_mc_nobck_unfold_effcor_cross_section");

  //  MnvH2D* disAMUMnv=(MnvH2D*)f2.Get("h_pzptrec_cross_section_qelike");
  //  MnvH2D* disNCTEQMnv=(MnvH2D*)f3.Get("h_pzptrec_cross_section_qelike");
  //  MnvH2D* ElasticFSIMnv=(MnvH2D*)f4.Get("h_pzptrec_cross_section_qelike");
  //  MnvH2D* MKMnv=(MnvH2D*)f5.Get("h_pzptrec_cross_section_qelike");
  MnvH2D* RPA_Res_MINOSMnv=(MnvH2D*)f6.Get("h_pzptrec_cross_section_qelike");
  //  MnvH2D* RPA_Res_NievesMnv=(MnvH2D*)f7.Get("h_pzptrec_cross_section_qelike");
  //  MnvH2D* ZExpansionMnv=(MnvH2D*)f8.Get("h_pzptrec_cross_section_qelike");
  MnvH2D* minerva_joint_lowq2Mnv=(MnvH2D*)f9.Get("h_pzptrec_cross_section_qelike");
    //  MnvH2D* defaultMnv=(MnvH2D*)f10.Get("h_pzptrec_cross_section_qelike");
  //  MnvH2D* pion_2p2hMnv=(MnvH2D*)f11.Get("h_pzptrec_cross_section_qelike");
  MnvH2D* pion_rpaMnv=(MnvH2D*)f12.Get("h_pzptrec_cross_section_qelike");
  //  MnvH2D* piontuneMnv=(MnvH2D*)f13.Get("h_pzptrec_cross_section_qelike");

  
  //External in TH3 format
  //NEUT
  TH3D *h_neut_540 = (TH3D*)f14.Get("plot");
  //  TH3D *h_neut_533_norpa = (TH3D*)f15.Get("plot");
  //  TH3D *h_neut_533_rpa = (TH3D*)f16.Get("plot");
  TH3D *h_neut_540_LFGNieves = (TH3D*)f17.Get("plot");

  //GENIE
  TH3D *h_genie_02a = (TH3D*)f18.Get("plot");
  TH3D *h_genie_02b = (TH3D*)f19.Get("plot");
  TH3D *h_genie_10a = (TH3D*)f20.Get("plot");
  TH3D *h_genie_10b = (TH3D*)f21.Get("plot");

  //NuWro
  TH3D *h_nuwro_sf = (TH3D*)f22.Get("plot");
  TH3D *h_nuwro_lfg = (TH3D*)f23.Get("plot");
  //GiBUU
  /*
  TH3D *h_gibuu = (TH3D*)f22.Get("ErecoilTIMES1E3VSmupzVSmuonpt_all");
  TH3D *h_gibuu_wc = (TH3D*)f23.Get("ErecoilTIMES1E3VSmupzVSmuonpt_all");

  TH3D *h_gibuu_qe = (TH3D*)f22.Get("ErecoilTIMES1E3VSmupzVSmuonpt_qe");
  TH3D *h_gibuu_res = (TH3D*)f22.Get("ErecoilTIMES1E3VSmupzVSmuonpt_res");
  TH3D *h_gibuu_dis = (TH3D*)f22.Get("ErecoilTIMES1E3VSmupzVSmuonpt_dis");
  TH3D *h_gibuu_2p2h = (TH3D*)f22.Get("ErecoilTIMES1E3VSmupzVSmuonpt_2p2h");
  TH3D *h_gibuu_other = (TH3D*)f22.Get("ErecoilTIMES1E3VSmupzVSmuonpt_other");
  */
 
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



  std::vector<std::vector<double> > full3D;
  full3D.push_back(recoil3Dbins);
  full3D.push_back(pt3Dbins);
  full3D.push_back(pz3Dbins);
  
  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);
  
  MnvH2D *mh_neut_540 = mcMnv->Clone("neut_540");
  //  MnvH2D *mh_neut_533_norpa = mcMnv->Clone("neut_533_norpa");
  //  MnvH2D *mh_neut_533_rpa = mcMnv->Clone("neut_533_rpa");
  MnvH2D *mh_neut_540_LFGNieves = mcMnv->Clone("neut_540_LFGNieves");

  MnvH2D *mh_genie_02a = mcMnv->Clone("genie_02a");
  MnvH2D *mh_genie_02b = mcMnv->Clone("genie_02b");
  MnvH2D *mh_genie_10a = mcMnv->Clone("genie_10a");
  MnvH2D *mh_genie_10b = mcMnv->Clone("genie_10b");

  MnvH2D *mh_nuwro_sf = mcMnv->Clone("nuwro_sf");
  MnvH2D *mh_nuwro_lfg = mcMnv->Clone("nuwro_lfg");
  /*
  MnvH2D *mh_genie_xianguo = mcMnv->Clone("genie_xianguo");

  MnvH2D *mh_gibuu = mcMnv->Clone("gibuu");
  MnvH2D *mh_gibuu_wc = mcMnv->Clone("gibuu_wc");
  MnvH2D *mh_gibuu_qe = mcMnv->Clone("gibuu_qe");
  MnvH2D *mh_gibuu_res = mcMnv->Clone("gibuu_res");
  MnvH2D *mh_gibuu_dis = mcMnv->Clone("gibuu_dis");
  MnvH2D *mh_gibuu_2p2h = mcMnv->Clone("gibuu_2p2h");
  MnvH2D *mh_gibuu_other = mcMnv->Clone("gibuu_other");
  */
  fillExternal(mh_neut_540, h_neut_540, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  fillExternal(mh_neut_540_LFGNieves, h_neut_540_LFGNieves ,  my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());

  fillExternal(mh_genie_02a, h_genie_02a, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  fillExternal(mh_genie_02b, h_genie_02b, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  fillExternal(mh_genie_10a, h_genie_10a, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  fillExternal(mh_genie_10b, h_genie_10b, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  cout << "NuWro" << endl;
  fillExternal(mh_nuwro_lfg, h_nuwro_lfg, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  fillExternal(mh_nuwro_sf, h_nuwro_sf, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());


  /*
  cout << "Xianguo" << endl;
  fillExternal(mh_genie_xianguo,h_genie_xianguo, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  cout << "Gibuu" << endl;
  fillExternal(mh_gibuu, h_gibuu, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  fillExternal(mh_gibuu_wc, h_gibuu_wc, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());


  fillExternal(mh_gibuu_qe, h_gibuu_qe, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  fillExternal(mh_gibuu_res, h_gibuu_res, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  fillExternal(mh_gibuu_dis, h_gibuu_dis, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  fillExternal(mh_gibuu_2p2h, h_gibuu_2p2h, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  fillExternal(mh_gibuu_other, h_gibuu_other, my3d, pt3Dbins.size(), pz3Dbins.size(), recoil3Dbins.size());
  */
  
  cout << "Writing to file" << endl;
  //Save externals to file for chi2 calc
  TFile *extout = new TFile("ExternalModels.root","RECREATE");
  mh_neut_540->Write();
  //  mh_neut_533_norpa->Write();
  //  mh_neut_533_rpa->Write();
  mh_neut_540_LFGNieves->Write();
  
  mh_genie_02a->Write();
  mh_genie_02b->Write();
  mh_genie_10a->Write();
  mh_genie_10b->Write();

  mh_nuwro_lfg->Write();
  mh_nuwro_sf->Write();
  /*
  mh_genie_xianguo->Write();
  mh_gibuu->Write();
  mh_gibuu_wc->Write();
  extout->Close();
  */
  std::cout << "Starting up getting the projections" << std::endl;
  std::vector<TH2D*> dataresults = my3d->Get2DHistos(dataMnv,true);
  std::vector<TH2D*> dataresults_statonly = my3d->Get2DHistos(dataMnv,false);
  std::vector<TH2D*> mcresults = my3d->Get2DHistos(mcMnv,false);


  //  std::vector<TH2D*> disAMU = my3d->Get2DHistos(disAMUMnv,false); 
  //  std::vector<TH2D*> disNCTEQ = my3d->Get2DHistos(disNCTEQMnv,false);
  //  std::vector<TH2D*> ElasticFSI = my3d->Get2DHistos(ElasticFSIMnv,false);
  //  std::vector<TH2D*> MK = my3d->Get2DHistos(MKMnv,false); 
  std::vector<TH2D*> RPA_Res_MINOS = my3d->Get2DHistos(RPA_Res_MINOSMnv,false); 
  //  std::vector<TH2D*> RPA_Res_Nieves = my3d->Get2DHistos(RPA_Res_NievesMnv,false); 
  //  std::vector<TH2D*> ZExpansion = my3d->Get2DHistos(ZExpansionMnv,false); 
  std::vector<TH2D*> minerva_joint_lowq2 = my3d->Get2DHistos(minerva_joint_lowq2Mnv,false); 
  //  std::vector<TH2D*> defaultresult = my3d->Get2DHistos(defaultMnv,false); 
  //  std::vector<TH2D*> pion_2p2h = my3d->Get2DHistos(pion_2p2hMnv,false); 
  std::vector<TH2D*> pion_rpa = my3d->Get2DHistos(pion_rpaMnv,false); 
  //  std::vector<TH2D*> piontune = my3d->Get2DHistos(piontuneMnv,false); 
  
  //External
  std::vector<TH2D*> neut_540 = my3d->Get2DHistos(mh_neut_540,false);
  //std::vector<TH2D*> neut_533_norpa = my3d->Get2DHistos(mh_neut_533_norpa,false);
  //std::vector<TH2D*> neut_533_rpa = my3d->Get2DHistos(mh_neut_533_rpa,false);
  std::vector<TH2D*> neut_540_LFGNieves = my3d->Get2DHistos(mh_neut_540_LFGNieves,false);

  std::vector<TH2D*> genie_02a = my3d->Get2DHistos(mh_genie_02a,false);
  std::vector<TH2D*> genie_02b = my3d->Get2DHistos(mh_genie_02b,false);
  std::vector<TH2D*> genie_10a = my3d->Get2DHistos(mh_genie_10a,false);
  std::vector<TH2D*> genie_10b = my3d->Get2DHistos(mh_genie_10b,false);


  std::vector<TH2D*> nuwro_lfg = my3d->Get2DHistos(mh_nuwro_lfg,false);
  std::vector<TH2D*> nuwro_sf = my3d->Get2DHistos(mh_nuwro_sf,false);
  /*
  std::vector<TH2D*> gibuu = my3d->Get2DHistos(mh_gibuu,false);
  std::vector<TH2D*> gibuu_wc = my3d->Get2DHistos(mh_gibuu_wc,false);

  std::vector<TH2D*> gibuu_qe = my3d->Get2DHistos(mh_gibuu_qe,false);
  std::vector<TH2D*> gibuu_res = my3d->Get2DHistos(mh_gibuu_res,false);
  std::vector<TH2D*> gibuu_dis = my3d->Get2DHistos(mh_gibuu_dis,false);
  std::vector<TH2D*> gibuu_2p2h = my3d->Get2DHistos(mh_gibuu_2p2h,false);
  std::vector<TH2D*> gibuu_other = my3d->Get2DHistos(mh_gibuu_other,false);
  */

  cout << "I found " << mcresults.size() << " histograms to work with " << endl;
  for(unsigned int i=1;i<pz3Dbins.size();i++){
    double pzwidth = pz3Dbins[i]-pz3Dbins[i-1];
    dataresults[i]->Scale(1e39/pzwidth,"width");
    dataresults_statonly[i]->Scale(1e39/pzwidth,"width");
    mcresults[i]->Scale(1e39/pzwidth,"width");

    //    disAMU[i]->Scale(1e39/pzwidth,"width");
    //    disNCTEQ[i]->Scale(1e39/pzwidth,"width");
    //    ElasticFSI[i]->Scale(1e39/pzwidth,"width");
    //    MK[i]->Scale(1e39/pzwidth,"width");
    RPA_Res_MINOS[i]->Scale(1e39/pzwidth,"width");
    //    RPA_Res_Nieves[i]->Scale(1e39/pzwidth,"width");
    //    ZExpansion[i]->Scale(1e39/pzwidth,"width");
    minerva_joint_lowq2[i]->Scale(1e39/pzwidth,"width");
    //    defaultresult[i]->Scale(1e39/pzwidth,"width");
    //    pion_2p2h[i]->Scale(1e39/pzwidth,"width");
    pion_rpa[i]->Scale(1e39/pzwidth,"width");
    //    piontune[i]->Scale(1e39/pzwidth,"width");

    
    //External
    neut_540[i]->Scale(1e39/pzwidth,"width");
    //    neut_533_norpa[i]->Scale(1e39/pzwidth,"width");
    //    neut_533_rpa[i]->Scale(1e39/pzwidth,"width");
    neut_540_LFGNieves[i]->Scale(1e39/pzwidth,"width");

    genie_02a[i]->Scale(1e39/pzwidth,"width");
    genie_02b[i]->Scale(1e39/pzwidth,"width");
    genie_10a[i]->Scale(1e39/pzwidth,"width");
    genie_10b[i]->Scale(1e39/pzwidth,"width");

    nuwro_sf[i]->Scale(1e39/pzwidth,"width");
    nuwro_lfg[i]->Scale(1e39/pzwidth,"width");
    /*
    genie_xianguo[i]->Scale(1e39/pzwidth,"width");
    gibuu[i]->Scale(1e39/(13.0*pzwidth),"width");
    gibuu_wc[i]->Scale(1e39/(13.0*pzwidth),"width");

    gibuu_qe[i]->Scale(1e39/(13.0*pzwidth),"width");
    gibuu_res[i]->Scale(1e39/(13.0*pzwidth),"width");
    gibuu_dis[i]->Scale(1e39/(13.0*pzwidth),"width");
    gibuu_2p2h[i]->Scale(1e39/(13.0*pzwidth),"width");
    gibuu_other[i]->Scale(1e39/(13.0*pzwidth),"width");
    */
  }

  vector<int> mycolors = getColors(2);
  vector<int> mycolors2 = getColors(4);
  for(unsigned int i=1;i<pz3Dbins.size();i++){//skip underflow and overflow pz bins
    cout << "DOING BIN " << i << endl;

    if(doRatio){
      dataresults[i]->Divide(mcresults[i]);
      dataresults_statonly[i]->Divide(mcresults[i]);

      //      disAMU[i]->Divide(mcresults[i]);
      //      disNCTEQ[i]->Divide(mcresults[i]);
      //      ElasticFSI[i]->Divide(mcresults[i]);
      //      MK[i]->Divide(mcresults[i]);
      RPA_Res_MINOS[i]->Divide(mcresults[i]);
      //      RPA_Res_Nieves[i]->Divide(mcresults[i]);
      //      ZExpansion[i]->Divide(mcresults[i]);
      minerva_joint_lowq2[i]->Divide(mcresults[i]);
      //      defaultresult[i]->Divide(mcresults[i]);
      //      pion_2p2h[i]->Divide(mcresults[i]);
      pion_rpa[i]->Divide(mcresults[i]);
      //      piontune[i]->Divide(mcresults[i]);
      
      //external
      neut_540[i]->Divide(mcresults[i]);
      //      neut_533_norpa[i]->Divide(mcresults[i]);
      //      neut_533_rpa[i]->Divide(mcresults[i]);
      neut_540_LFGNieves[i]->Divide(mcresults[i]);

      genie_02a[i]->Divide(mcresults[i]);
      genie_02b[i]->Divide(mcresults[i]);
      genie_10a[i]->Divide(mcresults[i]);
      genie_10b[i]->Divide(mcresults[i]);


      nuwro_lfg[i]->Divide(mcresults[i]);
      nuwro_sf[i]->Divide(mcresults[i]);
      /*
      genie_xianguo[i]->Divide(mcresults[i]);
      gibuu[i]->Divide(mcresults[i]);
      gibuu_wc[i]->Divide(mcresults[i]);

      gibuu_qe[i]->Divide(mcresults[i]);
      gibuu_res[i]->Divide(mcresults[i]);
      gibuu_dis[i]->Divide(mcresults[i]);
      gibuu_2p2h[i]->Divide(mcresults[i]);
      gibuu_other[i]->Divide(mcresults[i]);
      */
      //final
      mcresults[i]->Divide(mcresults[i]);

      
    }


      mcresults[i]->SetLineColor(kRed);
      mcresults[i]->SetLineWidth(2);

      // These line and marker styles will be propagated to the 1D plots
      dataresults[i]->SetMarkerStyle(kFullCircle);
      dataresults[i]->SetMarkerSize(0.7);
      dataresults[i]->SetLineColor(kBlack);
      dataresults[i]->SetLineWidth(2);

      dataresults_statonly[i]->SetMarkerStyle(1);
      dataresults_statonly[i]->SetMarkerSize(100);
      dataresults_statonly[i]->SetLineColor(kViolet);
      dataresults_statonly[i]->SetLineWidth(10);

      //      disAMU[i]->SetLineColor(mycolors[2]);
      //      disNCTEQ[i]->SetLineColor(mycolors[3]);
      //      ElasticFSI[i]->SetLineColor(mycolors[4]);
      //      MK[i]->SetLineColor(mycolors[5]);
      RPA_Res_MINOS[i]->SetLineColor(mycolors[6]);
      //      RPA_Res_Nieves[i]->SetLineColor(mycolors[7]);
      minerva_joint_lowq2[i]->SetLineColor(mycolors[8]);

      //      ZExpansion[i]->SetLineColor(mycolors[9]);
      //      defaultresult[i]->SetLineColor(mycolors[10]);
      //      pion_2p2h[i]->SetLineColor(mycolors[11]);
      pion_rpa[i]->SetLineColor(mycolors[12]);
      //      piontune[i]->SetLineColor(mycolors[13]);
      /*
      disAMU[i]->SetLineWidth(2);
      disNCTEQ[i]->SetLineWidth(2);
      ElasticFSI[i]->SetLineWidth(2);
      MK[i]->SetLineWidth(2);
      RPA_Res_MINOS[i]->SetLineWidth(2);
      RPA_Res_Nieves[i]->SetLineWidth(2);
      ZExpansion[i]->SetLineWidth(2);
      minerva_joint_lowq2[i]->SetLineWidth(2);
      defaultresult[i]->SetLineWidth(2);
      pion_2p2h[i]->SetLineWidth(2);
      pion_rpa[i]->SetLineWidth(2);
      piontune[i]->SetLineWidth(2);
      */
      //External
      neut_540[i]->SetLineColor(mycolors[13]);
      neut_540[i]->SetLineWidth(2);

      //      neut_533_norpa[i]->SetLineColor(mycolors[14]);
      //      neut_533_norpa[i]->SetLineWidth(2);

      //      neut_533_rpa[i]->SetLineColor(mycolors[15]);
      //      neut_533_rpa[i]->SetLineWidth(2);

      neut_540_LFGNieves[i]->SetLineColor(mycolors[16]);
      neut_540_LFGNieves[i]->SetLineWidth(2);

      
      genie_02a[i]->SetLineColor(mycolors[1]);
      genie_02a[i]->SetLineColor(mycolors[2]);
      genie_02a[i]->SetLineColor(mycolors[3]);
      genie_02a[i]->SetLineColor(mycolors[4]);

      genie_02a[i]->SetLineWidth(2);
      genie_02b[i]->SetLineWidth(2);
      genie_10a[i]->SetLineWidth(2);
      genie_10b[i]->SetLineWidth(2);

      nuwro_sf[i]->SetLineColor(mycolors2[1]);
      nuwro_sf[i]->SetLineWidth(2);
      nuwro_lfg[i]->SetLineColor(mycolors2[2]);
      nuwro_lfg[i]->SetLineWidth(2);
      /*
      gibuu[i]->SetLineColor(mycolors2[3]);
      gibuu[i]->SetLineWidth(2);

      gibuu_wc[i]->SetLineColor(mycolors2[4]);
      gibuu_wc[i]->SetLineWidth(2);


      gibuu_qe[i]->SetLineColor(mycolors[3]);
      gibuu_res[i]->SetLineColor(mycolors[4]);
      gibuu_dis[i]->SetLineColor(mycolors[5]);
      gibuu_2p2h[i]->SetLineColor(mycolors[6]);
      gibuu_other[i]->SetLineColor(mycolors[7]);
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
    histAndOpts.push_back(std::make_pair(dataresults[i], "histpe"));
    histAndOpts.push_back(std::make_pair(mcresults[i],       "graphlx"));

    if(set==0){
      //      histAndOpts.push_back(std::make_pair(defaultresult[i],       "histl"));
    }

    if(set==1){
      //      histAndOpts.push_back(std::make_pair(disAMU[i],       "histl"));
      //      histAndOpts.push_back(std::make_pair(disNCTEQ[i],       "histl"));
    }
    if(set==2){
      //      histAndOpts.push_back(std::make_pair(ElasticFSI[i],       "histl"));
    }
    if(set==3){
      //      histAndOpts.push_back(std::make_pair(MK[i],       "histl"));
    }
    if(set==4){
      histAndOpts.push_back(std::make_pair(RPA_Res_MINOS[i],       "histl"));
      //      histAndOpts.push_back(std::make_pair(RPA_Res_Nieves[i],       "histl"));
      histAndOpts.push_back(std::make_pair(minerva_joint_lowq2[i],       "histl"));
    }
    if(set==5){
      //      histAndOpts.push_back(std::make_pair(ZExpansion[i],       "histl"));
      //      histAndOpts.push_back(std::make_pair(pion_2p2h[i],       "histl"));
      histAndOpts.push_back(std::make_pair(pion_rpa[i],       "histl"));
      //      histAndOpts.push_back(std::make_pair(piontune[i],       "histl"));
    }
    
    if(set==6){
      
      histAndOpts.push_back(std::make_pair(nuwro_lfg[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(nuwro_sf[i],       "graphlx"));
      
    }
    if(set==7){
      histAndOpts.push_back(std::make_pair(genie_02a[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(genie_10a[i],       "graphlx"));
    }
    if(set==8){
      histAndOpts.push_back(std::make_pair(genie_02b[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(genie_10b[i],       "graphlx"));
    }
    if(set==9){
      histAndOpts.push_back(std::make_pair(genie_02a[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(genie_02b[i],       "graphlx"));
    }
    if(set==10){
      histAndOpts.push_back(std::make_pair(genie_10a[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(genie_10b[i],       "graphlx"));
    }

    if(set==11){
      histAndOpts.push_back(std::make_pair(neut_540[i],       "graphlx"));
      //      histAndOpts.push_back(std::make_pair(neut_533_norpa[i],       "graphlx"));
      //      histAndOpts.push_back(std::make_pair(neut_533_rpa[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(neut_540_LFGNieves[i],       "graphlx"));
    }
    /*
    if(set==18){
      histAndOpts.push_back(std::make_pair(gibuu[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(gibuu_wc[i],       "graphlx"));
    }
    if(set==19){
      histAndOpts.push_back(std::make_pair(gibuu[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(gibuu_qe[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(gibuu_res[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(gibuu_dis[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(gibuu_2p2h[i],       "graphlx"));
      histAndOpts.push_back(std::make_pair(gibuu_other[i],       "graphlx"));
    }
    */
    //    histAndOpts.push_back(std::make_pair(dataresults_statonly[i], "histpe"));
    //    histAndOpts.push_back(std::make_pair(dataresults[i],     "graph0 ep"));

    cout << "Mults" << endl;


    // ----------------------------------------------------------------------------------
    //
    // First make pt in bins of pz

    // Values to multiply each bin by to get them on a similar range
    //These are for the XY!! 
    vector<double> multipliers1 = doRatio?GetScalesRatio(histAndOpts,true,2.5):GetScales(histAndOpts,true,5.49,0.75,false);
    GridCanvas* gc=plotXAxis1D(histAndOpts, "#Sigma T_{p}","p_{t,#mu} (GeV)",  doMultipliers ? &multipliers1[0] : NULL);
    // Set the y range manually. Can also use gc->Remax() to guess automatically
    if(doRatio) gc->SetYLimits(0.01,2.01);
    else gc->SetYLimits(0.01, 5.49);
    if(doZoom) gc->SetXLimits(0,0.4);
    //gc->SetYLimits(0, 1e-39);
    gc->SetYTitle("d^{3}#sigma/dp_{T}dp_{||}d_{E_{vis}} (x10^{-39} cm^{2}/GeV^{3}/c^{3}/Nucleon)");
    //    gc->SetLogy(true);
    gc->Modified();
    // Example of adding a legend. The co-ordinate system is NDC on the
    // entire canvas, ie (0,0) in the bottom left corner of the canvas
    // (not the individual pad), and (1,1) in the top right
    TLegend* leg=new TLegend(0.8, 0.1, 1, 0.3);
    TLegend* bigleg = new TLegend(0,0,1,1);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(dataresults[i], "MINERvA data", "lpe");
    leg->AddEntry(mcresults[i], "MnvGENIE", "l");
    if(set==0){
			      //      leg->AddEntry(defaultresult[i], "MINERvA default GENIE", "l");
    }
    if(set==1){
      //      leg->AddEntry(disAMU[i], "MnvGENIE, disAMU", "l");
      //      leg->AddEntry(disNCTEQ[i], "MnvGENIE, disNCTEQ", "l");
    }
    if(set==2){
      //      leg->AddEntry(ElasticFSI[i], "MnvGENIE, ElasticFSI", "l");
    }
    if(set==3){
      //      leg->AddEntry(MK[i], "MnvGENIE, MK", "l");
    }
    if(set==4){
      leg->AddEntry(RPA_Res_MINOS[i], "MnvGENIE, RPA Res MINOS", "l");
      //      leg->AddEntry(RPA_Res_Nieves[i], "MnvGENIE, RPA Res Nieves", "l");
      leg->AddEntry(minerva_joint_lowq2[i], "MnvGENIE, minervalowq2", "l");
    }
    if(set==5){
      //      leg->AddEntry(ZExpansion[i], "MnvGENIE, ZExpansion", "l");
      //      leg->AddEntry(pion_2p2h[i], "Piontune+low recoil fit", "l");
      leg->AddEntry(pion_rpa[i], "Piontune+RPA", "l");
      //      leg->AddEntry(piontune[i], "Piontune", "l");
    }

    if(set==6){
      leg->AddEntry(nuwro_lfg[i], "NuWro LFG","l");
      leg->AddEntry(nuwro_sf[i], "NuWro SF","l");
    }

    if(set==7){
      leg->AddEntry(genie_02a[i], "GENIE 02a","l");
      leg->AddEntry(genie_10a[i], "GENIE 10a","l");
    }
    if(set==8){
      leg->AddEntry(genie_02b[i], "GENIE 02b","l");
      leg->AddEntry(genie_10b[i], "GENIE 10b","l");
    }
    if(set==9){
      leg->AddEntry(genie_02a[i], "GENIE 02a","l");
      leg->AddEntry(genie_02b[i], "GENIE 02b","l");
    }
    if(set==10){
      leg->AddEntry(genie_10a[i], "GENIE 10a","l");
      leg->AddEntry(genie_10b[i], "GENIE 10b","l");
    }
    if(set==11){
      leg->AddEntry(neut_540[i], "Neut 5.4.0,SF","l");
      //      leg->AddEntry(neut_533_norpa[i], "Neut 5.3.3 no RPA","l");
      //      leg->AddEntry(neut_533_rpa[i], "Neut 5.3.3 w/ RPA","l");
      leg->AddEntry(neut_540_LFGNieves[i], "Neut 5.4.0,LFG","l");
    }
    /*
    if(set==18){
      leg->AddEntry(gibuu[i], "GiBUU","l");
      leg->AddEntry(gibuu_wc[i], "GiBUU W&C","l");
    }
    if(set==19){
      leg->AddEntry(gibuu[i], "GiBUU","l");
      leg->AddEntry(gibuu_qe[i], "GiBUU-QE","l");
      leg->AddEntry(gibuu_res[i], "GiBUU-RES","l");
      leg->AddEntry(gibuu_dis[i], "GiBUU-DIS","l");
      leg->AddEntry(gibuu_2p2h[i], "GiBUU-2p2h","l");
      leg->AddEntry(gibuu_other[i], "GiBUU-Other","l");
    }
    */
    //set0
    bigleg->AddEntry(dataresults[i], "MINERvA data", "lpe");
    bigleg->AddEntry(mcresults[i], "MnvGENIE", "l");
    //    bigleg->AddEntry(defaultresult[i], "MINERvA default GENIE", "l");
    //set1
    //    bigleg->AddEntry(disAMU[i], "MnvGENIE, disAMU", "l");
    //    bigleg->AddEntry(disNCTEQ[i], "MnvGENIE, disNCTEQ", "l");
    //set2
    //    bigleg->AddEntry(ElasticFSI[i], "MnvGENIE, ElasticFSI", "l");
    //set3
    //    bigleg->AddEntry(MK[i], "MnvGENIE, MK", "l");
    //set4
    bigleg->AddEntry(RPA_Res_MINOS[i], "MnvGENIE, RPA Res MINOS", "l");
    //    bigleg->AddEntry(RPA_Res_Nieves[i], "MnvGENIE, RPA Res Nieves", "l");
    bigleg->AddEntry(minerva_joint_lowq2[i], "MnvGENIE, minervalowq2", "l");
    //set5
    //    bigleg->AddEntry(ZExpansion[i], "MnvGENIE, ZExpansion", "l");
    //    bigleg->AddEntry(pion_2p2h[i], "Piontune+low recoil fit", "l");
    bigleg->AddEntry(pion_rpa[i], "Piontune+RPA", "l");
    //    bigleg->AddEntry(piontune[i], "Piontune", "l");
    //set6
    /*
    bigleg->AddEntry(nuwro_lfg[i], "NuWro LFG","l");
    bigleg->AddEntry(nuwro_sf[i], "NuWro SF","l");
			    //set7

    bigleg->AddEntry(genie_ha_NievesQE_EmpMEC_NievesOther[i], "ha:nieves:emp", "l");
    bigleg->AddEntry(genie_ha_NievesQE_NievesMEC_NievesOther[i], "ha:nieves:nieves", "l");
    bigleg->AddEntry(genie_ha_NievesQE_SuSAMEC_NievesOther[i], "ha:nieves:susa", "l");
    
    bigleg->AddEntry(genie_ha_NievesQE_EmpMEC_NievesOtherFIXED[i], "ha:nieves:emp:FIXED", "l");
    bigleg->AddEntry(genie_ha_NievesQE_NievesMEC_NievesOtherFIXED[i], "ha:nieves:emp:FIXED", "l");
    bigleg->AddEntry(genie_ha_NievesQE_SuSAMEC_NievesOtherFIXED[i], "ha:nieves:emp:FIXED", "l");
    
    bigleg->AddEntry(genie_ha_LSQE_EmpMEC_NievesOther[i], "ha:lsqe:emp", "l");
    bigleg->AddEntry(genie_ha_LSQE_NievesMEC_NievesOther[i], "ha:lsqe:nieves", "l");
    bigleg->AddEntry(genie_ha_LSQE_SuSAMEC_NievesOther[i], "ha:lsqe:susa", "l");
    
    bigleg->AddEntry(genie_ha_LSQE_EmpMEC_NievesOtherFIXED[i], "ha:lsqe:emp:FIXED", "l");
    bigleg->AddEntry(genie_ha_LSQE_NievesMEC_NievesOtherFIXED[i], "ha:lsqe:nieves:FIXED", "l");
    bigleg->AddEntry(genie_ha_LSQE_SuSAMEC_NievesOtherFIXED[i], "ha:lsqe:susa:FIXED", "l");
    
    bigleg->AddEntry(genie_ha_SuSAQE_EmpMEC_NievesOther[i], "ha:susa:emp", "l");
    bigleg->AddEntry(genie_ha_SuSAQE_NievesMEC_NievesOther[i], "ha:susa:nieves", "l");
    bigleg->AddEntry(genie_ha_SuSAQE_SuSAMEC_NievesOther[i], "ha:susa:susa", "l");
    
    bigleg->AddEntry(genie_ha_SuSAQE_EmpMEC_NievesOtherFIXED[i], "ha:susa:emp:fixed", "l");
    bigleg->AddEntry(genie_ha_SuSAQE_NievesMEC_NievesOtherFIXED[i], "ha:susa:nieves:fixed", "l");
    bigleg->AddEntry(genie_ha_SuSAQE_SuSAMEC_NievesOtherFIXED[i], "ha:susa:susa:FIXED", "l");
    
    bigleg->AddEntry(genie_hn_NievesQE_SuSAMEC_NievesOther[i], "hn:nieves:susa", "l");
    bigleg->AddEntry(genie_hn_NievesQE_NievesMEC_NievesOther[i], "hn:nieves:nieves", "l");
    
    bigleg->AddEntry(genie_hn_NievesQE_SuSAMEC_NievesOtherFIXED[i], "hn:nieves:susa:fixed", "l");
    bigleg->AddEntry(genie_hn_NievesQE_NievesMEC_NievesOtherFIXED[i], "hn:nieves:nieves:fixed", "l");
    
    bigleg->AddEntry(genie_hn_SuSAQE_SuSAMEC_NievesOther[i], "hn:susa:susa", "l");
    bigleg->AddEntry(genie_hn_SuSAQE_NievesMEC_NievesOther[i], "hn:susa:nieves", "l");
    
    bigleg->AddEntry(genie_hn_SuSAQE_SuSAMEC_NievesOtherFIXED[i], "hn:susa:susa:fixed", "l");
    bigleg->AddEntry(genie_hn_SuSAQE_NievesMEC_NievesOtherFIXED[i], "hn:susa:nieves:fixed", "l");

    bigleg->AddEntry(genie_xianguo[i], "GENIE Xianguo Config","l");
    //set17
    bigleg->AddEntry(neut_540[i], "Neut 5.4.0,SF","l");
    bigleg->AddEntry(neut_533_norpa[i], "Neut 5.3.3 no RPA,LFG","l");
    bigleg->AddEntry(neut_533_rpa[i], "Neut 5.3.3 w/ RPA,LFG","l");
    bigleg->AddEntry(neut_540_LFGNieves[i], "Neut 5.4.0,LFG","l");
    //set18
    bigleg->AddEntry(gibuu[i], "GiBUU","l");
    bigleg->AddEntry(gibuu_wc[i], "GiBUU W&C","l");
    */

    TLatex mytex;
    mytex.SetTextSize(0.05);
    string mystring =     Form("%.2f < P_{||} [GeV] < %.2f",pz3Dbins[i-1],pz3Dbins[i]);
    mytex.DrawLatex(0.35,0.96,mystring.c_str());
    if(!doRatio)leg->Draw("SAME");
    if(doZoom){
      gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-pz-multiplier_bin_%d_set_%d-zoom.eps",i,set) : Form("nu-3d-xsec-comps-pz_bin_%d_set_%d-zoom.eps",i,set));
      gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-pz-multiplier_bin_%d_set_%d-zoom.png",i,set) : Form("nu-3d-xsec-comps-pz_bin_%d_set_%d-zoom.png",i,set));
      gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-pz-multiplier_bin_%d_set_%d-zoom.C",i,set) : Form("nu-3d-xsec-comps-pz_bin_%d_set_%d-zoom.C",i,set));
    }
    else{
      gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-pz-multiplier_bin_%d_set_%d.eps",i,set) : Form("nu-3d-xsec-comps-pz_bin_%d_set_%d.eps",i,set));
      gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-pz-multiplier_bin_%d_set_%d.png",i,set) : Form("nu-3d-xsec-comps-pz_bin_%d_set_%d.png",i,set));
      gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-pz-multiplier_bin_%d_set_%d.C",i,set) : Form("nu-3d-xsec-comps-pz_bin_%d_set_%d.C",i,set));
    }
    

    TCanvas *cleg = new TCanvas("Legend","Legend",10,10,1000,700);
    bigleg->Draw();
    cleg->Print("XSec_ModelComp_Legend.png");
    TCanvas *cleg2 = new TCanvas("Legend2","Legend2",10,10,1000,700);
    leg->Draw();
    cleg2->Print(Form("XSec_ModelComp_Legend_set_%d.png",set));
    cleg2->Print(Form("XSec_ModelComp_Legend_set_%d.eps",set));

  }
    
}

int main(int argc, char* argv[])
{

  string s_ratio = argv[2];
  string s_zoom = argv[3];
  bool doRatio = s_ratio=="1" ? true:false;
  bool doZoom = s_zoom=="1"?true:false;
  //linear
  for(int i=0;i<12;i++){
    //    if(i!=6) continue;
    if(!doRatio)makePlots(true,argv[1],false,doRatio,i,doZoom);
    else makePlots(false,argv[1],false,doRatio,i,doZoom);

  //log
    //    if(!doRatio) makePlots(true,argv[1],true,doRatio,i,doZoom);
    //    else makePlots(false,argv[1],true,doRatio,i,doZoom);
  }

  return 0;
}
