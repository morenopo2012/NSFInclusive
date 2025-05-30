#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
//#include "syst_common.h"
#include "TArrayD.h"
#include <iostream>
                              
#include "Math/Factory.h"     
#include "Math/Functor.h"     
#include "Math/Minimizer.h"   

#include "RooUnfold/RooUnfoldBayes.h"
#include "RooUnfold/RooUnfoldSvd.h"
#include "RooUnfold/RooUnfoldTUnfold.h"
#include "RooUnfold/RooUnfoldInvert.h"
#include "RooUnfold/RooUnfoldBinByBin.h"
#include "RooUnfold/RooUnfoldResponse.h"
#include "RooUnfold/RooUnfold.h"
#include "MinervaUnfold/MnvResponse.h"
#include "MinervaUnfold/MnvUnfold.h"

#include "PlotUtils/MacroUtil.h" 
//#include "include/CVUniverse_faiza.h"
#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "../include/Variable.h"
//#include "../../include/Variable_faiza.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "../../NUKECCSRC/ana_common/include/LateralSystematics.h"
#include <iostream>
#include <stdlib.h>
//#include "Cintex/Cintex.h"
#include "TParameter.h"
#include "../include/plotting_functions.h"
#ifndef __CINT__
#endif

using namespace NUKECC_ANA;
using namespace PlotUtils;
using namespace std;

namespace {bool isoCorrect = false;
          
  double ratio_min = .5;
  double ratio_max = 2.5;
           }

void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name);
void GetMinAndMaxAxisRatio( TH1D* histData, TH1D* histMC, double dataMCScale, double &plotMin, double &plotMax );



//const bool do_fits = true;//Run this code with this boolean as true first and then run it again with boolean as false before plotiing
const bool do_fits = false;//second

TH1D* GetVertErrorBandUniverseHist(MnvH1D* mnvh1d, const std::string& errName, int universe); 

double getPhysChi2( const double * par);
      TH1D* m_histo_trans_data;
      TH1D* m_histo_sig_trans;
      TH1D* m_histo_trans_trans;
      TH1D* m_histo_contin_trans;
      TH1D* m_histo_contin_data;
      TH1D* m_histo_sig_contin;
      TH1D* m_histo_trans_contin;
      TH1D* m_histo_contin_contin;


 
int UnfoldIterations( const std::string& var );


int main(int argc, char * argv[]){
  //ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if(argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n"<<
       "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
      "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" <<
      "\t-Playlist\t \t =\t eg. minervame1A"<< std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    return 0;
  }

  string outdir=argv[1];
  int targetID = atoi(argv[2]);
  int targetZ = atoi(argv[3]);
  const string playlist= argv[4];

  TString hname, hname1;
   if(RunCodeWithSystematics){
  hname1= Form("%s/Hists_EventSelection_t%d_z%02d_Nu.root", outdir.c_str(), targetID, targetZ);
  hname= Form("%s/Unfolded_t%d_z%02d_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str());
     }
   else{
  hname= Form("%s/Hists_PhysicsBackgd_without_SYS_t%d_z%02d_Nu_v1.root", outdir.c_str(), targetID, targetZ);
     }
  TFile *f2 = new TFile( hname1,"read" );
  TFile *f1 = new TFile( hname,"read" );
  cout<<hname<<endl;

  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(playlist);
  NukeCC_Binning  *binsDef = new NukeCC_Binning();
 
  MnvPlotter mnvPlotter(kNukeCCStyle);

  string trueZ, mat;
  if( targetZ == 26 ){
      trueZ = "Iron";
      mat = "Fe";
   }
  if( targetZ == 82 ){
       trueZ = "Lead";
       mat = "Pb";
   }
  if( targetZ == 6 ){
        trueZ = "Carbon";
        mat = "C";
   }

  const string suffix = Form( "t%02d_z%02d_%s", targetID, targetZ, trueZ.c_str() );

  vector<string> vars;
  vars.clear();
  vars.push_back("Emu");
  vars.push_back("Enu");
  vars.push_back("Ehad");
  //vars.push_back("x");
  //vars.push_back("y");

   // cout<<"Am I here?????"<<endl;
   
   HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(playlist);
   //TString unfoldFileName = utils->GetHistFileName( "Unfolded", FileType::kAny, targetID, targetZ, helicity );
  /* TFile *fUnfold = new TFile( unfoldFileName, "RECREATE" );
   assert(fUnfold);

   TString fileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ, helicity );
   TFile *fin = new TFile( fileName, "RECREATE" );
  

   TString BGfileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ, helicity );
   TFile *BGfin = new TFile( BGfileName, "RECREATE" );
*/
///////////////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/////////////

////  TFile *fRawXSec = new TFile( utils->GetXSecFileName().Data() );
   
    //cout<<"Am I here?????"<<endl;

  TFile *fxsec = new TFile( Form("%s/CrossSection_z%d_%s.root", outdir.c_str(), targetZ, playlist.c_str() ), "recreate" );
  //TFile *fxsec= new TFile(fxsecname, "RECREATE");
    assert(fxsec);    
//  TString fxsecname= utils->GetHistFileName("CrossSection", FileType::kAny, 0, targetZ, helicity);
//  TFile *fxsec= new TFile(fxsecname, "RECREATE");

    
  double prelimX = 0.7, prelimY = 0.7;

///////////////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/////////////

  TParameter<double> *mcPOT = (TParameter<double>*)f2->Get("MCPOT");
  TParameter<double> *dataPOT = (TParameter<double>*)f2->Get("DataPOT");
  double mcpot = mcPOT->GetVal();
  double datapot = dataPOT->GetVal();

 
  double dataMCScale = datapot/mcpot;
  cout<<"Data POT = "<<datapot<<"   MC POT = "<<mcpot<<"   Data/MC POT = "<<dataMCScale<<endl;

    // get the hists
    std::vector<MnvH1D*> hists, hists_acccorr, hists_dis;

    vector< MnvH1D*> hists_trans_data, hists_contin_data;
    vector< MnvH1D* > hists_true_signal_in_trans,  hists_true_signal_in_contin;
    vector< MnvH1D* > hists_true_trans_in_trans,   hists_true_trans_in_contin;
    vector< MnvH1D* > hists_true_contin_in_contin, hists_true_contin_in_trans;
    vector< MnvH1D* > hists_true_signal_in_signal;
    vector< MnvH1D* > hists_trans_mc_CH_US,  hists_trans_mc_CH_DS;
    vector< MnvH1D* > hists_contin_mc_CH_US, hists_contin_mc_CH_DS;
    vector< MnvH1D* > tmp_scale_factor_trans, tmp_scale_factor_contin;
 
    vector< MnvH1D* >true_trans_in_trans_tuned;
    vector< MnvH1D* >true_contin_in_trans_tuned;
    vector< MnvH1D* >true_trans_in_contin_tuned;
    vector< MnvH1D* >true_contin_in_contin_tuned;
    vector< MnvH1D* >tuned_trans_mc_US, tuned_trans_mc_DS, tuned_contin_mc_US, tuned_contin_mc_DS; 

    vector<MnvH1D*> hists_true_wrong_target_USCH, hists_true_wrong_target_DSCH, hists_true_wrong_target_Other;
    vector<MnvH1D*> hists_true_dis_wrong_target_USCH, hists_true_dis_wrong_target_DSCH, hists_true_dis_wrong_target_Other;
    std::vector<MnvH1D*> Hists, mcSignals, mcHistsBG, mcHists_nondis, FullSysts;
    std::vector<MnvH1D*> hists_true_sig, hists_trans_SB, hists_contin_SB;
  
  //Get the cross section hists
  vector<MnvH1D*> GENIExsecHists, GENIEEvtHists;
  cout<<"Getting the cross section histogram from "<<utils->GetXSecFileName()<<endl;
  

  for( vector<string>::iterator i = vars.begin(); i != vars.end(); ++i )
    {
       const string& var = *i;

      cout<<"plotting variable  "<<var<<endl;
      f1->cd();      
      MnvH1D *dSum(0);
      MnvH1D *mcSum(0);
      //MnvH1D * dSum  = utils->GetXSecRecoHistAvg( var, FileType::kAny, targetZ, helicity, isoCorrect, f1  );
      //MnvH1D * mcSum = utils->GetXSecRecoHistAvg( var, FileType::kAny, targetZ, helicity, isoCorrect, f1 );
      
      fxsec->cd();
      
      ////dSum->Write(Form("effcorrected_data_%s", var.c_str() ));
      ////mcSum->Write(Form("effcorrected_mc_%s", var.c_str() ));
     cout<<"Am I here?????"<<endl; 
      const string defaultCut = "std_dis";
      MnvH1D *xsecPassive(0);
      //        MnvH1D *xsecPassive = NukeUtils::Get().GetXSecHistAvg( fRawXSec, targetZ, var, defaultCut, isoCorrect);
            
      if( !dSum )// || !xsecPassive)
	{
	  cout << " no data Could not get passive target sum for var = " << var << endl;
	  continue;
	}

      if( !mcSum )// || !xsecPassive)
	{
	  cout << " no mc Could not get passive target sum for var = " << var << endl;
	  continue;
	}
  
      for( int bin=0;bin<dSum->GetNbinsX();bin++){
	cout<<bin<< " data content  "<<dSum->GetBinContent(bin)<<"  "<<dSum->GetBinError(bin)<<"  "<<dSum->GetBinError(bin)/dSum->GetBinContent(bin) <<endl;
	cout<<bin<< " mc content    "<<mcSum->GetBinContent(bin)<<"  "<<mcSum->GetBinError(bin)<<"  "<<mcSum->GetBinError(bin)/mcSum->GetBinContent(bin) <<endl;
	//            cout<<bin<< " genie content "<<xsecPassive->GetBinContent(bin)<<endl;
      }
      
      MnvH1D *dSumXSec  = (MnvH1D*)dSum->Clone("tmpXSec_Passive_D");
      MnvH1D *mcSumXSec = (MnvH1D*)mcSum->Clone("tmpXSec_Passive_MC");
            
      MnvH1D *dSum_EffCorr = (MnvH1D*)dSum->Clone("tmp_EffCorr_Passive_D");
      MnvH1D *mcSum_EffCorr = (MnvH1D*)mcSum->Clone("tmp_EffCorr_Passive_MC");
           
      std::string tmp_playlist;
      if(playlist=="AllNuME") tmp_playlist = "minervame1A";
      else tmp_playlist=playlist;
   

      if( "Enu" == var )
	{
                
	  utils->DivideByFlux( dSumXSec,  tmp_playlist, true );
	  utils->DivideByFlux( mcSumXSec, tmp_playlist, true );
	}
      else
	{
	  utils->DivideByIntegratedFlux( dSumXSec, tmp_playlist, true );
	  utils->DivideByIntegratedFlux( mcSumXSec, tmp_playlist, true );
	}
            
      double scattering_centers_data = utils->GetTotalScatteringCenters( targetZ, 0 );
      double scattering_centers_mc   = utils->GetTotalScatteringCenters( targetZ, 1 );
            
      cout<<"Number of scattering centers MC "<<scattering_centers_mc<<" data "<<scattering_centers_data<<endl;


      if(targetZ==99){
	scattering_centers_data = scattering_centers_data*9.0;
	scattering_centers_mc   = scattering_centers_mc*9.0;
      }
            
      if(playlist=="AllNuME"){

	mcSumXSec->Scale( 1E38/ (datapot*scattering_centers_mc)   );
	dSumXSec->Scale(  1E38/ (datapot*scattering_centers_data) );
	// xsecPassive->Scale(  1E38);
      }
      else{
	mcSumXSec->Scale( 1E38/ (mcpot*scattering_centers_mc)     );
	dSumXSec->Scale(  1E38/ (datapot*scattering_centers_data) );
	// xsecPassive->Scale(  1E38 );
      }
            
      dSumXSec->Write(Form("xsec_data_%s", var.c_str() ));
      mcSumXSec->Write(Form("xsec_mc_%s", var.c_str() ));

      MnvH1D *mcXSecPerE =  (MnvH1D*)mcSumXSec->Clone("tmpXSecPerE_MC");
      mcXSecPerE->Reset();
      for(int bin = 1; bin < mcXSecPerE->GetNbinsX() + 1; ++bin){
	double binE     = mcSumXSec->GetBinCenter(bin);
	double xsection = mcSumXSec->GetBinContent(bin);
	mcXSecPerE->SetBinContent(bin, xsection/binE);
                
      }
      mcXSecPerE->Write(Form("xsecPerE_mc_%s", var.c_str() ));
	    
      const double minVal = binsDef->GetDISVarMinVal( var );
      const double maxVal = binsDef->GetDISVarMaxVal( var );
      dSum->GetXaxis()->SetRangeUser( minVal, maxVal );
      mcSum->GetXaxis()->SetRangeUser( minVal, maxVal );
            
      dSumXSec->GetXaxis()->SetRangeUser( minVal, maxVal );
      mcSumXSec->GetXaxis()->SetRangeUser( minVal, maxVal );
      mcXSecPerE->GetXaxis()->SetRangeUser( minVal, maxVal );


      if( NO_SPREAD_ERROR )
	{
	  dSum->SetUseSpreadErrorAll(false);
	  mcSum->SetUseSpreadErrorAll(false);
	  dSumXSec->SetUseSpreadErrorAll(false);
	  mcSumXSec->SetUseSpreadErrorAll(false);
	}
            
            
      if( xsecPassive )
	{
	  xsecPassive->GetXaxis()->SetRangeUser( minVal, maxVal );
	  xsecPassive->SetMarkerColor( kBlue - 2 );
	  xsecPassive->SetLineColor( kBlue - 2 );
	  xsecPassive->SetMarkerStyle( 34 );
	}

      string x_title = binsDef->GetXaxisTitle( var );
      if( utils->UnfoldVar( var ) )
	x_title = "Unfolded " + x_title;
      else
	x_title = "Reconstructed " + x_title;
      const string y_title = binsDef->GetYaxisTitle( var ) + " / nucleon";
      const string xsec_name  = var == "Enu" ?
        utils->GetTotalXSecString( targetID, targetZ ):
	utils->GetXSecString( var, targetID, targetZ );
      const string xsec_units = utils->GetTotalXSecUnits();
            
      double potX = .31, potY = .88;

            
      if(var=="Enu"){
	dSumXSec->SaveAs(Form("DataXSec_%d.root",targetZ));
	mcSumXSec->SaveAs(Form("MCXSec_%d.root",targetZ));
      }


/*
      //====== Sum of Passive ======/
      {
	const string y_label = "#sigma" + xsec_units;
                
	TString cName = Form( "SumResultPlot_CrossSectionDataMC_%s_%s", var.c_str(), suffix.c_str() );
	TCanvas c( cName, cName, 1280, 800 );
	mcSumXSec->GetXaxis()->SetTitle( x_title.c_str() );
	mcSumXSec->GetYaxis()->SetTitle( y_label.c_str() );
	//mnvPlotter.axis_maximum = 35.0;
	//            mnvPlotter.DrawDataMCWithErrorBand( dSumXSec, mcSumXSec, dataMCScale, "TC-TR");
	cout<<cName<<endl;
	//DrawDataMCWithErrorBand(         datahist,    mchist, mcscale, legpos, usehisttitles, bgdhist, databkghist, covareanormalize, statplussys,issmooth, mcstatonly)
	//This one is actual function. Has to modify PlotUtils
	//mnvPlotter.DrawDataMCWithErrorBand( dSumXSec, mcSumXSec, 1.0, "TC-TR", false, NULL, NULL, false, true, false, true);
	mnvPlotter.DrawDataMCWithErrorBand( dSumXSec, mcSumXSec, 1.0, "TC-TR", false, NULL, NULL, false);
	//At the moment this is drawing with full systmatic uncertainties on the MC which is wrong. Need to add some option in here. 
	if( WRITE_PRELIMINARY )
	  mnvPlotter.WritePreliminary(.32,.87);
	if( WRITE_POT ) mnvPlotter.WriteNorm("", potX, potY, .05, datapot );
	if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( xsec_name.c_str(), .04 );
	mnvPlotter.MultiPrint( &c );
                
                
                
	cName = Form( "SumResultPlot_CrossSectionDataMCRatio_%s_%s", var.c_str(), suffix.c_str() );
	TCanvas cR( cName, cName, 1280, 800 );
	mcSumXSec->GetYaxis()->SetTitle( "Data/MC" );
	//            mnvPlotter.DrawDataMCRatio( dSumXSec, mcSumXSec, dataMCScale, true, true, ratio_min, ratio_max );
	////mnvPlotter.DrawDataMCRatio( dSumXSec, mcSumXSec, 1.0, false, true, true, ratio_min, ratio_max );
	mnvPlotter.DrawDataMCRatio( dSumXSec, mcSumXSec, 1.0, true, true, ratio_min, ratio_max );
	//            mnvPlotter.AddChi2Label( dSumXSec, mcSumXSec, dataMCScale, "TC", .06 );
	mnvPlotter.AddChi2Label( dSumXSec, mcSumXSec, 1.0, "TC", .06 );
	if( WRITE_PRELIMINARY )
	  mnvPlotter.WritePreliminary(.7,.78);
	if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( xsec_name.c_str(), .04 );
	mnvPlotter.MultiPrint( &cR );
                
	cName = Form( "SumResultPlot_EffCorrectedDist_%s_%s",var.c_str(), suffix.c_str() );
	TCanvas cEff( cName, cName, 1280, 800);
	mcSum_EffCorr->GetXaxis()->SetTitle( x_title.c_str() );
	mcSum_EffCorr->GetYaxis()->SetTitle("Efficiency Corrected Events");
	//mnvPlotter.DrawDataMCWithErrorBand(dSum_EffCorr, mcSum_EffCorr, dataMCScale, "TC-TR",true);
	//DrawDataMCWithErrorBand(         datahist,    mchist, mcscale, legpos, usehisttitles, bgdhist, databkghist, covareanormalize, statplussys,issmooth, mcstatonly)
	//need to modify PlotUtils to incorporate this function
	//mnvPlotter.DrawDataMCWithErrorBand( dSum_EffCorr, mcSum_EffCorr, dataMCScale, "TC-TR", true, NULL, NULL, false, true, false, true);
	mnvPlotter.DrawDataMCWithErrorBand( dSum_EffCorr, mcSum_EffCorr, dataMCScale, "TC-TR", true, NULL, NULL, false, true);
	if( WRITE_PRELIMINARY )
	  mnvPlotter.WritePreliminary(.32,.87);
	if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( xsec_name.c_str(), .04 );
	mnvPlotter.MultiPrint( &cEff );
                
                
	cName = Form( "SumResultPlot_EffCorrectedDist_DataErrSummary_%s_%s",var.c_str(), suffix.c_str() );
	TCanvas cEffUnc( cName, cName, 1280, 800);
	dSum_EffCorr->GetXaxis()->SetTitle( x_title.c_str() );
	dSum_EffCorr->GetYaxis()->SetTitle("Fractional Uncertainty");
	mnvPlotter.DrawErrorSummary(dSum_EffCorr,"TL",true);
	if( WRITE_PRELIMINARY )
	  mnvPlotter.WritePreliminary(.32,.87);
	if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( xsec_name.c_str(), .04 );
	mnvPlotter.MultiPrint( &cEffUnc );
    }
//========== End of Passive Sum ============//
      if(xsecPassive)
	{
	  const int old_dataMarker = mnvPlotter.data_marker;
	  const int old_dataColor  = mnvPlotter.data_color;
	  const double old_dataSize = mnvPlotter.data_marker_size;
	  mnvPlotter.data_marker = 34;
	  mnvPlotter.data_color  = kBlue-2;
	  mnvPlotter.data_marker_size = 2.5;
                
	  const string y_label = "#sigma" + xsec_units;
	  TString cName = Form("SumResultPlot_CrossSectionMCGenerated_%s_%s", var.c_str(), suffix.c_str() );
	  TCanvas c( cName, cName, 1280, 800 );
	  mcSumXSec->SetTitle("MC");
	  xsecPassive->SetTitle("GENIE Calculation");
	  mcSumXSec->GetXaxis()->SetTitle( x_title.c_str() );
	  mcSumXSec->GetYaxis()->SetTitle( y_label.c_str() );
	  //            mnvPlotter.DrawDataMCWithErrorBand( xsecPassive, mcSumXSec, dataMCScale, "TC-TR", true);
	  mnvPlotter.DrawDataMCWithErrorBand( xsecPassive, mcSumXSec, 1.0, "TC-TR", true);
	  if( WRITE_PRELIMINARY )
	    mnvPlotter.WritePreliminary(.32,.87);
	  if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( xsec_name.c_str(), .04 );
	  mnvPlotter.MultiPrint( &c );
                
	  double oldHeadroom = mnvPlotter.headroom;
	  mnvPlotter.headroom = 1.2;
                
	  cName = Form( "SumResultPlot_CrossSectionMCGeneratedRatio_%s_%s", var.c_str(), suffix.c_str() );
	  TCanvas cR( cName, cName, 1280, 800 );
	  xsecPassive->GetXaxis()->SetTitle( x_title.c_str() );
	  xsecPassive->GetYaxis()->SetTitle( "MC / GENIE Calculation" );
	  xsecPassive->GetYaxis()->SetRangeUser( 0.4, 1.2 );
	  //            mnvPlotter.DrawDataMCRatio( mcSumXSec, xsecPassive, 1./dataMCScale, false, true, 0.15, 1.4, "MC / GENIE Calculation" );
	  //This is actual one I have to use in future
	  //mnvPlotter.DrawDataMCRatio( mcSumXSec, xsecPassive, 1., false, false, true, 0.15, 1.4, "MC / GENIE Calculation" ); //I might end up with error bands I don't want here now... cause I lost the ability to turn them off. 
	  mnvPlotter.DrawDataMCRatio( mcSumXSec, xsecPassive, 1., false, false, 0.15, 1.4 ); //I might end up with error bands I don't want here now... cause I lost the ability to turn them off. 
	  //            mnvPlotter.AddChi2Label( (TH1*)xsecPassive, (TH1*)mcSumXSec, dataMCScale, "TC", .06 );
	  mnvPlotter.AddChi2Label( (TH1*)xsecPassive, (TH1*)mcSumXSec, 1.0, "TC", .06 );
	  if( WRITE_PRELIMINARY )
	    mnvPlotter.WritePreliminary("TR");
	  if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( xsec_name.c_str(), .04 );
	  mnvPlotter.MultiPrint( &cR );
                
	  //return data marker style
	  mnvPlotter.data_marker = old_dataMarker;
	  mnvPlotter.data_color  = old_dataColor;
	  mnvPlotter.data_marker_size = old_dataSize;
	  mnvPlotter.headroom = oldHeadroom;
                
         
                
	  //            MnvH1D *dSum_EffCorr = (MnvH1D*)dSum->Clone("tmp_EffCorr_Passive_D");
	  //            MnvH1D *mcSum_EffCorr = (MnvH1D*)mcSum->Clone("tmp_EffCorr_Passive_MC");
                
                
	}
          
*/
            
      //===============================
      // ERRORS
      //===============================
            
      mnvPlotter.height_nspaces_per_hist = 1.;
      mnvPlotter.legend_n_columns = 2;
      mnvPlotter.headroom=1.6;
            
      string legLoc = "TL";
      string prelimLoc = "TR";
            
      //=======================
      // Passive sum
      //=======================
      {
	TString cName = Form( "SumResultPlot_CrossSectionErrSummary_Data_%s_%s", var.c_str(), suffix.c_str() );
	TCanvas cErrSummary( cName, cName, 1280,920 );
	mcSum->GetXaxis()->SetTitle( x_title.c_str() );
	mcSum->GetYaxis()->SetTitle( "Fractional Uncertainty" );
	bool allSolid = true;
	//            mnvPlotter.DrawDataMCErrorSummary( mcSum, dSum, legLoc, false, allSolid );
	//             mnvPlotter.DrawDataMCErrorSummary( dSumXSec, mcSumXSec, legLoc, false, allSolid );
	mnvPlotter.DrawErrorSummary(dSumXSec,legLoc,true);
                
	if( WRITE_PRELIMINARY )
	  mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
	if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Uncertainties on %s", xsec_name.c_str() ), .05);
	mnvPlotter.MultiPrint( &cErrSummary );
      }
            
      {
	TString cName = Form( "SumResultPlot_CrossSectionErrSummaryMC_%s_%s", var.c_str(), suffix.c_str() );
	TCanvas cErrSummary( cName, cName, 1280,920 );
	mcSum->GetXaxis()->SetTitle( x_title.c_str() );
	mcSum->GetYaxis()->SetTitle( "Fractional Uncertainty" );
	bool allSolid = true;
	//            mnvPlotter.DrawDataMCErrorSummary( mcSum, dSum, legLoc, false, allSolid , 0.0001, false, "MC Stats." );
	//            mnvPlotter.DrawDataMCErrorSummary( dSumXSec, mcSumXSec, legLoc, false, allSolid , 0.0001, false );
	mnvPlotter.DrawErrorSummary(mcSumXSec,legLoc,true);
	if( WRITE_PRELIMINARY )
	  mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
	if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Uncertainties on %s", xsec_name.c_str() ), .05);
	mnvPlotter.MultiPrint( &cErrSummary );
      }
            
      {
	TString cName = Form( "SumResultPlot_CrossSectionErrSummary_Data_DetRes_%s_%s", var.c_str(), suffix.c_str() );
	TCanvas cErrSummary( cName, cName, 1280,920 );
	mcSum->GetXaxis()->SetTitle( x_title.c_str() );
	mcSum->GetYaxis()->SetTitle( "Fractional Uncertainty" );
	bool allSolid = true;
	mnvPlotter.DrawErrorSummary( dSumXSec, "TC", false, true, 0.0001, false, "Detector Res." );
	//mnvPlotter.DrawDataMCErrorSummary( dSumXSec, mcSumXSec, legLoc, false, allSolid , 0.0001, false, "Detector Res." );
	if( WRITE_PRELIMINARY )
	  mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
	if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Uncertainties on %s", xsec_name.c_str() ), .05);
	mnvPlotter.MultiPrint( &cErrSummary );
      }
            
      {
	TString cName = Form( "SumResultPlot_CrossSectionErrSummary_Data_FSIModels_%s_%s", var.c_str(), suffix.c_str() );
	TCanvas cErrSummary( cName, cName, 1280,920 );
	mcSum->GetXaxis()->SetTitle( x_title.c_str() );
	mcSum->GetYaxis()->SetTitle( "Fractional Uncertainty" );
	bool allSolid = true;
	mnvPlotter.DrawErrorSummary( dSumXSec, "TC", false, true, 0.0001, false, "FSI Models" );
	//            mnvPlotter.DrawDataMCErrorSummary( dSumXSec, mcSumXSec, legLoc, false, allSolid , 0.0001, false, "FSI Models" );
	if( WRITE_PRELIMINARY )
	  mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
	if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Uncertainties on %s", xsec_name.c_str() ), .05);
	mnvPlotter.MultiPrint( &cErrSummary );
      }
            
            
      mnvPlotter.height_nspaces_per_hist = 1.2;
      mnvPlotter.legend_n_columns = 1;
            
      //cleanup
      delete mcSum;
      delete dSum;
      if( xsecPassive )
	delete xsecPassive;
            
      mnvPlotter.CleanTmp();
    }
        
//  delete fRawXSec;
  delete fxsec;        
  return 0;


  
    

}//variable loop closed








double getPhysChi2( const double * par )
{
   
    double scale_trans  = par[0];
    double scale_contin = par[1];
    
    // Sideband is bins from plane number 11 to 65
    // Assume no error on the MC, and use the stat errors on the data
    int binLow  = m_histo_trans_data->FindBin( 2.0  ); // same bins on both sets of histograms
    int binHigh = m_histo_trans_data->FindBin( 25.0 );
    
    double chi2 = 0.0;
    // clone and scale the floating histograms
    TH1D * temp_trans_float = new TH1D( *(TH1D*)m_histo_trans_trans->Clone("temp_trans_float") );
    TH1D * temp_contin_trans_float = new TH1D( *(TH1D*)m_histo_contin_trans->Clone("temp_contin_trans_float") );
    
    temp_trans_float->Scale( scale_trans );
    temp_contin_trans_float->Scale( scale_contin );
    temp_trans_float->Add(temp_contin_trans_float);
    temp_trans_float->Add(m_histo_sig_trans);
    
    TH1D * temp_contin_float = new TH1D( *(TH1D*)m_histo_contin_contin->Clone("temp_contin_float") );
    TH1D * temp_trans_contin_float = new TH1D( *(TH1D*)m_histo_trans_contin->Clone("temp_trans_contin_float") );
    temp_contin_float->Scale( scale_contin );
    temp_trans_contin_float->Scale( scale_trans );
    temp_contin_float->Add(temp_trans_contin_float);
    temp_contin_float->Add(m_histo_sig_contin);
    
    double tot_trans_data   = 0.0;
    double tot_trans_float  = 0.0;
    double tot_contin_data  = 0.0;
    double tot_contin_float = 0.0;
    
    for( int b = binLow; b <= binHigh; ++b ) {
        double trans_data   = m_histo_trans_data->GetBinContent(b);
        double trans_float  = temp_trans_float->GetBinContent(b);
        
        if( trans_data != 0.0 ){
            tot_trans_data  += trans_data;
            tot_trans_float += trans_float;
            
            chi2 += pow( (trans_float - trans_data)/sqrt(trans_data),2);
        }
    }
    
    for( int b = binLow; b <= binHigh; ++b ) {
        double contin_data   = m_histo_contin_data->GetBinContent(b);
        double contin_float  = temp_contin_float->GetBinContent(b);
        
        if( contin_data != 0.0 ){
            tot_contin_data  += contin_data;
            tot_contin_float += contin_float;
            
            chi2 += pow( (contin_float - contin_data)/sqrt(contin_data),2);
        }
    }
    
    double unbinned_chi2 = pow( (tot_contin_float-tot_contin_data)/sqrt(tot_contin_data), 2 );
    unbinned_chi2 += pow( (tot_trans_float-tot_trans_data)/sqrt(tot_trans_data), 2 );
    
    
//    cout<<"The chi2 integrated is "<<unbinned_chi2<<"  and the unintegrated is "<<chi2<<endl;
    return unbinned_chi2;
    
     //return chi2;
    
}





void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name ){
  
  //first declared some plotting parameters
  mnvPlotter.legend_offset_x = .03;
  mnvPlotter.width_xspace_per_letter = .33;
  /*
  string targetString;  
  if( targetZ == 6 )
    targetString = "Carbon";
  else if( targetZ == 26 )
    targetString = "Iron";
  if( targetZ == 82 )
    targetString = "Lead";
 */
   
  //just get the CV for now 
  MnvH1D* data_cv = (MnvH1D*)dataHisto->Clone("mnvh1d_data_histo");
  cout << "data_cv = " << data_cv->Integral() << endl; 
  //get total MC
  MnvH1D* totalMC = (MnvH1D*)mchistos[0]->Clone("mnvh1d_mc_histo");
  totalMC->UseCurrentStyle();

  /*string cname = (string)canname;
  if( targetZ == 6 && string::npos != cname.find("_US_") && string::npos != cname.find("_planeDNN_") ) totalMC->GetXaxis()->SetRangeUser(22., 28.);
  else if( targetZ == 6 && string::npos != cname.find("_DS_") && string::npos != cname.find("_planeDNN_") ) totalMC->GetXaxis()->SetRangeUser(34., 40.);
  */
  for( int i = 1; i < mchistos.GetLast()+1; i++ ){
    MnvH1D* mc_histo_add = (MnvH1D*)mchistos[i]->Clone("mnvh1d_mc_histo_add");
    mc_histo_add->SetFillColor(0);
    mc_histo_add->SetFillStyle(0);
    
    totalMC->Add(mc_histo_add);
  }
  
  cout << "total MC = " << totalMC->Integral() << endl; 
  if( data_cv->GetSumw2N() == 0 ) 
    data_cv->Sumw2();
  if( totalMC->GetSumw2N() == 0 ){ 
    totalMC->Sumw2();
  }


  
  int ndfStat = 1;
  double chi2Sys = mnvPlotter.Chi2DataMC( data_cv, totalMC, ndfStat, dataMCScale );
  ndfStat -= 1;
  
  TString labelStat = Form("Stat. + Sys. Error, #chi^{2}/ndf = %3.2f/%d = %3.2f", chi2Sys, ndfStat, chi2Sys/(Double_t)ndfStat);
  
 // string region = "Fiducial";
 // if( string::npos != cname.find("_US_") ) region = "Upstream Plastic Background";
 // if( string::npos != cname.find("_DS_") ) region = "Downstream Plastic Background";
  
 // string status = "After Tuning";
 // if( string::npos != cname.find("BeforeTuning_") ) status = "Before Tuning";

 // if( writeChi2 && !tune_fiducial ) fprintf(Chi2Table, "%s %s %s %.2f \n", var.c_str(), region.c_str(), status.c_str(), chi2Sys/(Double_t)ndfStat ); 
 // if( writeChi2 && tune_fiducial  ) fprintf(Chi2Table, "%s Fiducial %s, With Sys. %.2f \n", var.c_str(), status.c_str(), chi2Sys/(Double_t)ndfStat ); 

  double plotMin = -1., plotMax = -1.;
  GetMinAndMaxAxisRatio( data_cv, totalMC, dataMCScale, plotMin, plotMax );
  
  //plot data/mc with stat errors
  {
  TString cName;
  if(plotUS) cName = Form("Ratio_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
  if(plotDS) cName = Form("Ratio_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
    //TString cName = Form("BGPlotRatio_%s",var.c_str());
    //TString Yaxis = "Data/MC (" + status + ")"; 
    TString Yaxis = "Data/MC"; 
    std::string x_label = var.c_str();
    TCanvas c( cName, cName, 1200, 800 );
    cout << "plotting " << cName << endl;
    cout << "data_cv =" << data_cv->GetBinContent(2) << endl;
    cout << "mc_cv =" << totalMC->GetBinContent(2) << endl;
    //mnvPlotter.DrawDataMCRatio( data_cv, totalMC, dataMCScale, true, true, plotMin, plotMax, Yaxis );
    mnvPlotter.DrawDataMCRatio( data_cv, totalMC, dataMCScale, true, true, 0., 2., x_label.c_str(), Yaxis );
    mnvPlotter.AddPlotLabel( labelStat, .46, .875, .05 );
    cout<<"Is this correct number = "<<dataMCScale<<endl;
 
    if( WRITE_PRELIMINARY ) mnvPlotter.WritePreliminary("BC");
    
    //if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("%s of %s", region.c_str(), targetString.c_str() ), .06);
    mnvPlotter.MultiPrint( &c );
  }
  
  {    
    TString cName;
   if(plotUS) cName = Form("FracError_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
   if(plotDS) cName = Form("FracError_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
    TCanvas can( cName, cName, 1200,800 );
    
    cout << "Fill color: "<< totalMC->GetFillColor() << endl;
 
    can.SetFillColor(0);
    can.SetFillStyle(4050);
    totalMC->UseCurrentStyle();
    
  string legLoc  = "TL";
  string prelimLoc = "TR";
  std::string x_label = var.c_str(); 
    cout << "plotting " << cName << endl;
    //totalMC->PopVertErrorBand(ErrorName::PlasticBG);
    //cout<<"Canname:" <<canname<<"\t totalMC Plastic_SB? "<< totalMC->HasUncorrError(ErrorName::CHSB)<<endl;
    //mnvPlotter.DrawDataMCErrorSummary( totalMC, data_cv, /*cname,*/ legLoc, false, true, 0.00001, false ); //use my colors and styles
    PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
    //PlotUtils::MnvPlotter plotter(kCCQENuStyle);
    //plotter.DrawErrorSummary(totalMC, legLoc);
    //mnvPlotter.DrawErrorSummary( totalMC, legLoc, false, true, 0.0001, true, "", "", "", false); 
    plotter.DrawErrorSummary( totalMC, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
    
    
    if(WRITE_PRELIMINARY)
      mnvPlotter.WritePreliminary(prelimLoc);
    
    mnvPlotter.MultiPrint( &can, can.GetName() );
  }
  //scale to data
  //totalMC->Scale(dataMCScale);
}

void GetMinAndMaxAxisRatio( TH1D* histData, TH1D* histMC, double dataMCScale, double &plotMin, double &plotMax ){
  
  TH1D* ratio = (TH1D*)histData->Clone("ratio_plot");
  TH1D* totalMC = (TH1D*)histMC->Clone("mc_clone");
  totalMC->Scale(dataMCScale);
  ratio->Divide( histData, totalMC);
  
  int firstBin = ratio->GetXaxis()->GetFirst();
  int lastBin = ratio->GetXaxis()->GetLast();
  
  std::map<double, int> bins;
  for( int b = firstBin; b <= lastBin; b++ ){
    bins.insert( std::pair<double,int>(ratio->GetBinContent(b), b) ); 
  }
  
  vector<int> bins_in_order;
  for( std::map<double,int>::iterator it = bins.begin(); it != bins.end(); ++it ){
    bins_in_order.push_back( it->second );
  }
  
  int last_bin = bins_in_order.size() - 1;
  plotMax = ratio->GetBinContent( bins_in_order[last_bin] ) + ratio->GetBinError( bins_in_order[last_bin] ) + 0.3;
  if( plotMax > 4.0 ) plotMax = ratio->GetBinContent( bins_in_order[last_bin-1] ) + ratio->GetBinError( bins_in_order[last_bin-1] ) + 0.1;
  plotMin = ratio->GetBinContent( bins_in_order[0] ) - ratio->GetBinError( bins_in_order[0] ) - 0.1;
  if( plotMin <= -0.1 ) plotMin = ratio->GetBinContent( bins_in_order[1] ) - ratio->GetBinError( bins_in_order[1] ) - 0.1;
  
}


//MnvH1D* ReBinPhysSBTuning( MnvH1D* h_scaleFactors, string SB, int targetID, int targetZ, string playlist, string var, bool cv_only, NukeCC_Binning *binsDef,  MnvH1D* h_sb){
MnvH1D* ReBinPhysSBTuning( MnvH1D* h_scaleFactors, string SB, int targetID, int targetZ, string playlist, string var, bool cv_only,  MnvH1D* h_sb){
    cout<<"I'm saving the other binnings for the sidebands now "<<endl;
    int temp_tgtZ = 0;
    if(targetID > 10)
        temp_tgtZ = 99;
    if(targetID< 10)
        temp_tgtZ = targetZ;
    NukeCC_Binning  *binsDef = new NukeCC_Binning();     
    int nonzerobin = h_scaleFactors->FindFirstBinAbove(0.);

    vector<double> binsLowEdge = binsDef->GetEnergyBins( var, true ); //Do I really want this? The sideband histograms are using SidebandBins (only relevant for Q2 and W) probably?
    vector<double> SBbins = binsDef->GetSidebandBins( var ); 
    TString hName = Form("scale_factor_%s_%s", SB.c_str(),var.c_str());
    MnvH1D* tunedHisto = new MnvH1D( hName, hName, SBbins.size()-1, &SBbins[0] ); //binsLowEdge.size()-1, &binsLowEdge[0] ); 
    MnvH1D* scaleHisto = (MnvH1D*)tunedHisto->Clone(Form("scale_factor_%s_%s", SB.c_str(),var.c_str()));
    scaleHisto->Reset();
    scaleHisto->ClearAllErrorBands();
    //cout<<" I cloned the histogram "<<endl;
    double cv_val = h_scaleFactors->GetBinContent(nonzerobin);
    cout<<"I've got the scale factor stuff for SB "<<SB.c_str()<<"  the CV value is "<<cv_val<<endl;

    
    for(int i=0;i<scaleHisto->GetNbinsX()+2;i++) scaleHisto->SetBinContent(i,cv_val);
    
    std::vector<std::string> vertNames = h_scaleFactors->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k ){
        MnvVertErrorBand *errBand = h_scaleFactors->GetVertErrorBand( vertNames[k] );
        const int universes = errBand->GetNHists();
        std::vector<TH1D*> vert_hists;
        for(int u=0;u<universes;++u){
            TH1D* tmp_scale = new TH1D(*errBand->GetHist( u ));
            TH1D* tmp_template = new TH1D(scaleHisto->GetCVHistoWithStatError());
            tmp_template->SetName(Form("Scale_universe_%d",u));
            for(int i=0;i<scaleHisto->GetNbinsX()+2;i++) tmp_template->SetBinContent(i,tmp_scale->GetBinContent(nonzerobin));
            vert_hists.push_back(tmp_template);
        }
        scaleHisto->AddVertErrorBand( vertNames[k],vert_hists);
    }
    // Uncorrelated errors
    if(h_sb != NULL){
      int offset=0;
      if(binsLowEdge[0] != SBbins[0]){
	cout<<"The EnergyBins and SidebandBins differ."<<endl;
	/*for(int j=0; j<=h_sb->GetNbinsX()+2; ++j){
	  if (binsLowEdge[0] == SBbins[j]){
	    offset=j;
	    cout<<"Var: "<<var.c_str()<<"\t offset: "<<j<< "Value: "<<SBbins[j]<<endl;
	    break;
	  }
	}//end loop over bins
	if(offset==0){
            Error( "NukeUtils::ReBinPhysSBTuning", "Inconsistent binnings between EnergyBins and SidebandBins. Check that SidebandBins contains EnergyBins" );
            throw 1;
	}
	*/
      }//end if case for different binnings
    cout<<"What about here"<<endl;
      std::vector<std::string> uncorrNames = h_sb->GetUncorrErrorNames();
      scaleHisto->AddMissingErrorBandsAndFillWithCV(*h_sb);
      
      for (std::vector<std::string>::iterator name = uncorrNames.begin(); name != uncorrNames.end(); ++name){
	std::string errName = *name;
	TH1D* sferr= scaleHisto->GetUncorrError(errName);
	TH1D* tmp_err= h_sb->GetUncorrError(errName); 
	

	for( int i=0; i<=sferr->GetNbinsX()+2;++i){//keep the fractional error from the untuned sideband
	  sferr->SetBinContent(i,cv_val); //make sure CV is consistent
	  double mccontent=tmp_err->GetBinContent(i);
	  double mcerr=tmp_err->GetBinError(i);
	  if(mccontent==0.0) sferr->SetBinError(i, 0.0);
	  else sferr->SetBinError(i, cv_val * mcerr/mccontent);
	  
	} 

      }//end loop over uncorr errors
    }// end if h_sb not NULL



    cout<<"Scale histogram CV "<<scaleHisto->GetBinContent(2)<<endl;
    return scaleHisto;
//    scaleHisto->Write();
//    scaleFactor->Close();
}
