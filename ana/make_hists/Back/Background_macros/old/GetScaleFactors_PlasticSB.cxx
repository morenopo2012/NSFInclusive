#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "include/NukeCCUtilsNSF.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
//#include "syst_common.h"
#include "TArrayD.h"
#include <iostream>
                              
#include "Math/Factory.h"     
#include "Math/Functor.h"     
#include "Math/Minimizer.h"   


//#include "include/CVUniverse_faiza.h"
#include "include/CVUniverse.h"
#include "../../include/Variable.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "include/LateralSystematics.h"
#include <iostream>
#include <stdlib.h>

#ifndef __CINT__
#endif

using namespace NUKECC_ANA;
using namespace PlotUtils;
using namespace std;

  const bool makeDISPlots = false; //set to true to make DIS plots!
  const bool makeLowQ2Plots = false; //set to true to make DIS plots!
  const bool makeLowWPlots = false; //set to true to make DIS plots!
  const bool onlyCalculateCV = false; //set to false to do the tuning in systematics universe - set to true if you're not calculating systematics
  const bool writeSysCode = true; //set to true to write the systematics scale factors code

  const bool writeChi2 = false;
  const bool tune_fiducial = false; //set to true to see effects in passive targets
  FILE *Chi2Table;

  //plotting parameters
  double statX = .68, statY = .44;
  string legPos = "TL", normPos="TR"; 
  string legLoc  = "TL";
  string prelimLoc = "TR";
  double infoX = .68;
  double infoY = .78;
  

//----Helper Functions-----------------------------------------------
//void SetHistosAxis( MnvH1D* histo, string var, const bool makeDISPlots );
void SetStackedHistosMaxAxis( MnvPlotter &mnvPlotter, int targetZ, string var, const bool makeDISPlots, int opt );
void CalcScaleFactorAndTuneBG( MnvPlotter mnvPlotter, vector<MnvH1D*> histos_fit, vector<MnvH1D*> histos_mc_tuning, const string var, int targetID, int targetZ, double dataMCScale, TFile *scaleFactor, bool fixUS, bool fixDS, string playlist, bool makeDISPlots);
//! Return a string to add to the printed plot
//TString GetPlotSuffix() const;
//void TuneFiducial( MnvPlotter &mnvPlotter, vector<MnvH1D*> histos_mc_tuning, int targetZ, string var, TString canname, string plot_title, double dataMCScale, string playlist );
//void SubtractFiducial( MnvPlotter &mnvPlotter, vector<MnvH1D*> histos_mc_tuning, int targetZ, string var, TString canname, string plot_title, double dataMCScale, string playlist );
void ConstrainScintBackground(MnvH1D* mnvh1d_data, MnvH1D* mnvh1d_mc_signal, MnvH1D* mnvh1d_mc_other, MnvH1D* mnvh1d_mc_background_US, MnvH1D* mnvh1d_mc_background_DS,MnvH1D* mc_signal, MnvH1D* mc_other, MnvH1D* mc_background_US, MnvH1D* mc_background_DS, MnvH1D* tmp_scale_factor, bool fixUS, bool fixDS, string channel, bool onlyCV , bool writeScaleFactor );
vector<double> CalcScaleFactorMinimizer( TH1D* h1d_data, TH1D* h1d_mc_signal, TH1D* h1d_mc_other, TH1D* h1d_mc_backgroundUS, TH1D* h1d_mc_backgroundDS, const std::string& errName, bool fixUS, bool fixDS);
double getChi2( const double * par );
      TH1D* m_histo_data;
      TH1D* m_histo_sig_fix;
      TH1D* m_histo_other_fix;
      TH1D* m_histo_plastic_fix;
      TH1D* m_histo_plastic_float;
      TH1D* m_histo_plastic_us_float;
      TH1D* m_histo_plastic_ds_float;
void AddHistosToAVectorUS( vector<MnvH1D*> &v_histos, vector<MnvH1D*> histos, int targetZ ); 
void AddHistosToAVectorDS( vector<MnvH1D*> &v_histos, vector<MnvH1D*> histos, int targetZ ); 
TH1D* GetVertErrorBandUniverseHist(MnvH1D* mnvh1d, const std::string& errName, int universe); 
void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetZ, TString canname );
void DrawMeSomeDataMCStackedHistos( MnvPlotter mnvPlotter, MnvH1D* histo_data, TObjArray hlist, double scale, TString canname, string plotname, string var, int targetZ, const bool makeDISPlots );
//void DrawFracMCHistos( MnvPlotter mnvPlotter, TObjArray hlist, TString canname, string plotname, string var, int targetZ, const bool makeDISPlots );
//-------------------------------------------------------------------

int main(int argc, char * argv[]){
  ROOT::Cintex::Cintex::Enable();
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

   TString histFileName;
   if(RunCodeWithSystematics){
       histFileName = Form("%s/Hists_PlasticBackgd_with_SYS_t%d_z%02d_Nu_v1_NuMAD.root", outdir.c_str(), targetID, targetZ ); 
     }
   else{
       histFileName = Form("%s/Hists_PlasticBackgd_without_SYS_t%d_z%02d_Nu_v1_NuMAD.root", outdir.c_str(), targetID, targetZ);
     } 
 
  cout<<histFileName<<endl;

  vector< MnvH1D* > hists_sideband_US, hists_sideband_DS;

  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(playlist);

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

  vector<string> vars;
  vars.clear();
  vars.push_back("planeDNN");
  vars.push_back("Ehad");
  vars.push_back("Enu");
  vars.push_back("Emu");
//  vars.push_back("Q2");
//  vars.push_back("W");
//  vars.push_back("x");
//  vars.push_back("y");


  TFile *f1 = new TFile( histFileName,"read" );

  vector< MnvH1D*> hists_US_data, hists_US_mat, hists_US_other, hists_US_regUS, hists_US_regDS;
  vector< MnvH1D*> hists_DS_data, hists_DS_mat, hists_DS_other, hists_DS_regUS, hists_DS_regDS;
  vector< MnvH1D*> hists_tgt_data, hists_tgt_mat, hists_tgt_other, hists_tgt_regUS, hists_tgt_regDS;
  vector< MnvH1D*> histos_fit_US, histos_mc_tuning_US, histos_fit_DS, histos_mc_tuning_DS, histos_mc_tuning_fid;

 for( vector<string>::iterator i = vars.begin(); i != vars.end(); ++i )
    {
       const string& var = *i;

  f1->cd();
   hists_US_data.push_back(  (MnvH1D*)f1->Get(Form("selected_data_reco_US_%s", var.c_str())));
   hists_US_mat.push_back(   (MnvH1D*)f1->Get(Form("US_%s_%s_%s", mat.c_str(), trueZ.c_str(), var.c_str())));
   hists_US_other.push_back( (MnvH1D*)f1->Get(Form("US_other_%s_%s", trueZ.c_str(), var.c_str())));
   hists_US_regUS.push_back( (MnvH1D*)f1->Get(Form("US_regUS_%s_%s", trueZ.c_str(), var.c_str())));
   hists_US_regDS.push_back( (MnvH1D*)f1->Get(Form("US_regDS_%s_%s", trueZ.c_str(), var.c_str())));
  
   hists_DS_data.push_back(  (MnvH1D*)f1->Get(Form("selected_data_reco_DS_%s", var.c_str())));
   hists_DS_mat.push_back(   (MnvH1D*)f1->Get(Form("DS_%s_%s_%s", mat.c_str(), trueZ.c_str(), var.c_str())));
   hists_DS_other.push_back( (MnvH1D*)f1->Get(Form("DS_other_%s_%s", trueZ.c_str(), var.c_str())));
   hists_DS_regUS.push_back( (MnvH1D*)f1->Get(Form("DS_regUS_%s_%s", trueZ.c_str(), var.c_str())));
   hists_DS_regDS.push_back( (MnvH1D*)f1->Get(Form("DS_regDS_%s_%s", trueZ.c_str(), var.c_str())));
 
   hists_tgt_data.push_back(  (MnvH1D*)f1->Get(Form("selected_data_reco_tgt_%s", var.c_str())));
   hists_tgt_mat.push_back(   (MnvH1D*)f1->Get(Form("Tgt_%s_%s_%s",  mat.c_str(), trueZ.c_str(), var.c_str())));
   hists_tgt_other.push_back( (MnvH1D*)f1->Get(Form("Tgt_other_%s_%s", trueZ.c_str(), var.c_str())));
   hists_tgt_regUS.push_back( (MnvH1D*)f1->Get(Form("Tgt_regUS_%s_%s", trueZ.c_str(), var.c_str())));
   hists_tgt_regDS.push_back( (MnvH1D*)f1->Get(Form("Tgt_regDS_%s_%s", trueZ.c_str(), var.c_str())));
}

//  double datapot = 8.98e+19;//utils->getPOTData(f1);
//  double mcpot = 2.00e+20;//utils->getPOTMC(f1);
  double datapot = 1.06091e+21;//  whole Nu data set
  double mcpot = 2.51803e+21;//  whole Nu data set
 
  double dataMCScale = datapot/mcpot;
  cout<<"Data POT = "<<datapot<<"   MC POT = "<<mcpot<<"   Data/MC POT = "<<dataMCScale<<endl;
 
/* 
  double datapot = utils->getPOTData(f1);
  double mcpot = utils->getPOTMC(f1);
  double dataMCScale = datapot/mcpot;
  cout<<"Data POT = "<<datapot<<"   MC POT = "<<mcpot<<"   Data/MC POT = "<<dataMCScale<<endl;
*/ 
for( unsigned int i = 0; i != 1; ++i )
    {
       const string& var = vars[i];
       cout<<" variable in the planeDNN loop = "<<var<<endl;
 if(var == "planeDNN"){
  histos_fit_US.push_back( hists_US_data[i] );
  histos_fit_US.push_back( hists_US_mat[i]  );
  histos_fit_US.push_back( hists_US_other[i] );
  histos_fit_US.push_back( hists_US_regUS[i] );
  histos_fit_US.push_back( hists_US_regDS[i] );
  
  histos_fit_DS.push_back( hists_DS_data[i] );
  histos_fit_DS.push_back( hists_DS_mat[i] );
  histos_fit_DS.push_back( hists_DS_other[i] );
  histos_fit_DS.push_back( hists_DS_regUS[i] );
  histos_fit_DS.push_back( hists_DS_regDS[i] );
    }
  else{
   cout<<"not taking zeorth compement"<<endl;
    }
  }

  TFile *scaleFactor;
      
  if( ! tune_fiducial ){
    scaleFactor = new TFile( Form("%s/CHScaleFactors_t%d_z%d_%s_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str(), getenv("NUKECC_TAG") ), "recreate" );
    assert(scaleFactor);
  }
 

for( unsigned int i = 0; i != vars.size(); ++i )
    {
       const string& var = vars[i];
       cout<<"******************************"<<endl;
       cout<<"       variable = "<<var<<endl;
       cout<<"******************************"<<endl;
     
  histos_mc_tuning_US.push_back( hists_US_data[i] );
  histos_mc_tuning_US.push_back( hists_US_mat[i] );
  histos_mc_tuning_US.push_back( hists_US_other[i] );
  histos_mc_tuning_US.push_back( hists_US_regUS[i] );
  histos_mc_tuning_US.push_back( hists_US_regDS[i] );
      
  histos_mc_tuning_DS.push_back( hists_DS_data[i] );
  histos_mc_tuning_DS.push_back( hists_DS_mat[i] );
  histos_mc_tuning_DS.push_back( hists_DS_other[i] );
  histos_mc_tuning_DS.push_back( hists_DS_regUS[i] );
  histos_mc_tuning_DS.push_back( hists_DS_regDS[i] );
  
   scaleFactor->cd();
      
    //pass the vectors of histos for background fitting (USING MINUIT2MINIMIZER) and tuning
    CalcScaleFactorAndTuneBG( mnvPlotter, histos_fit_US, histos_mc_tuning_US, vars[i], targetID, targetZ, dataMCScale, scaleFactor, false, true, playlist, makeDISPlots );
      
    CalcScaleFactorAndTuneBG( mnvPlotter, histos_fit_DS, histos_mc_tuning_DS, vars[i], targetID, targetZ, dataMCScale, scaleFactor, true, false, playlist, makeDISPlots );
      
    //clear vector before going to the next variable 
    histos_mc_tuning_US.clear();
    histos_mc_tuning_DS.clear();
      
    /*    } else {
	    
      //All Selected DIS - Downstream Region
      histos_mc_tuning_fid.push_back( hists_data_fiducial[ivar] );
      AddHistosToAVector( histos_mc_tuning_fid, hists_fiducial[ivar], targetZ );
      
      //pass the vectors of histos for background fitting and tuning
      string plot_title = "Tuned Fiducial -  " + trueZ;
      canname = Form( "%s_%s",  vars[ivar].c_str(), suffix.c_str() );
      TuneFiducial( mnvPlotter, histos_mc_tuning_fid, targetZ, vars[ivar], canname, plot_title, dataMCScale, playlist );
      
      //pass the vectors of histos for background fitting and tuning
      plot_title = "BG Subtracted -  " + trueZ;
      canname = Form( "%s_%s",  vars[ivar].c_str(), suffix.c_str() );
      SubtractFiducial( mnvPlotter, histos_mc_tuning_fid, targetZ, vars[ivar], canname, plot_title, dataMCScale, playlist );
      
      //clear vector before going to the next variable
      histos_mc_tuning_fid.clear(); //second index is DS scale
      }
      */
        }//loop over all variables
      
      //write scale factor histograms 
      if( ! tune_fiducial ) scaleFactor->Write();
      if( ! tune_fiducial ) scaleFactor->Close();
      
      //close text file containing chi^2/ndf values
      if( writeChi2 )fclose(Chi2Table);
      
      //*********************************************************END TUNING STUFF}
      return 0;
}


void  CalcScaleFactorAndTuneBG( MnvPlotter mnvPlotter, vector<MnvH1D*> histos_fit, vector<MnvH1D*> histos_mc_tuning, const string var, int targetID, int targetZ, double dataMCScale, TFile *scaleFactor, bool fixUS, bool fixDS, string playlist, bool makeDISPlots){
  
  //only want to save scale factors histograms
  string channelstr = "Inclusive";
  if(makeDISPlots) channelstr = "DIS";
 
  //clone histograms for fitting
  MnvH1D* data_histo        = (MnvH1D*)histos_fit[0]->Clone("histo_data");
  MnvH1D* signal_histo     = (MnvH1D*)histos_fit[1]->Clone("histo_mc_signal");
  MnvH1D* other_histo       = (MnvH1D*)histos_fit[2]->Clone("histo_mc_other_targ");
  MnvH1D* plastic_us_histo  = (MnvH1D*)histos_fit[3]->Clone("histo_mc_plasticUS");
  MnvH1D* plastic_ds_histo  = (MnvH1D*)histos_fit[4]->Clone("histo_mc_plasticDS"); 
 
  //scale Data/MC
    signal_histo->Scale(dataMCScale);
    other_histo->Scale(dataMCScale);
    plastic_us_histo->Scale(dataMCScale);
    plastic_ds_histo->Scale(dataMCScale);
          
  //clone for before tuned mc templates
  MnvH1D* data_sideband_untuned       = (MnvH1D*)histos_mc_tuning[0]->Clone("untuned_data");
  MnvH1D* untuned_mc_sideband_signal  = (MnvH1D*)histos_mc_tuning[1]->Clone("untuned_mc_signal");
  MnvH1D* untuned_mc_sideband_Other   = (MnvH1D*)histos_mc_tuning[2]->Clone("untuned_mc_other");
  MnvH1D* untuned_mc_sideband_regUS      = (MnvH1D*)histos_mc_tuning[3]->Clone("untuned_mc_us");
  MnvH1D* untuned_mc_sideband_regDS      = (MnvH1D*)histos_mc_tuning[3]->Clone("untuned_mc_us");
  //clone the tuned mc templates
  MnvH1D* data_sideband             = (MnvH1D*)histos_mc_tuning[0]->Clone("tuned_data");
  MnvH1D* tuned_mc_sideband_signal = (MnvH1D*)histos_mc_tuning[1]->Clone("tuned_mc_signal");
  MnvH1D* tuned_mc_sideband_Other   = (MnvH1D*)histos_mc_tuning[2]->Clone("tuned_mc_other");
  MnvH1D* tuned_mc_sideband_regUS      = (MnvH1D*)histos_mc_tuning[3]->Clone("tuned_mc_us");
  MnvH1D* tuned_mc_sideband_regDS      = (MnvH1D*)histos_mc_tuning[4]->Clone("tuned_mc_ds");
 
  //untunned scintillator background before adding it to the stacks 
  TObjArray hlistUntuned;
  hlistUntuned.Add(untuned_mc_sideband_signal);
  hlistUntuned.Add(untuned_mc_sideband_Other);
  hlistUntuned.Add(untuned_mc_sideband_regUS);
  hlistUntuned.Add(untuned_mc_sideband_regDS);
  
  //---- PLOT DATA/MC RATIO and CHI^2
  string targetString;
  if( targetZ == 6 ) targetString = "t03_z06";
  if( targetZ == 26 ) targetString = "t01_z26";
  if( targetZ == 82 ) targetString = "t01_z82";
  
/*
  //Figure out what this was supposed to be doing before I deleted the canvases
  string cname = (string)canname;
  cout <<"CNAME = "<<cname<<endl;// Defining name for the plot  "AfterTuning_....".
  string regID;
  if( string::npos != cname.find("_US_") ){
    regID = "US";
  }
  if( string::npos != cname.find("_DS_") ){
    regID = "DS";
    }
    
  TString cName = "BeforeTuning_" + regID + "_" + var + "_" + targetString; // + "_" + GetPlotSuffix().Data(); 
  cout<<"cName = "<<cName<<endl;
  PlotChi2Stat( mnvPlotter, data_sideband_untuned, hlistUntuned, dataMCScale, var, targetZ, cName ); 
*/ 
  MnvH1D* tmp_scale_factor;
  TString histName;
  if( !fixUS ){ // if fixUS is false => fixUS = 0, i.e. it will take US events. 
    histName = "scaleFactor_US_" + var; 
    tmp_scale_factor = (MnvH1D*)histos_mc_tuning[3]->Clone(histName);
  }
  if( !fixDS ){ // if fixDS is false => fixDS = 0, i.e. it will take DS events. 
    histName = "scaleFactor_DS_" + var; 
    tmp_scale_factor = (MnvH1D*)histos_mc_tuning[4]->Clone(histName);
  }

  ConstrainScintBackground( data_histo, signal_histo, other_histo, plastic_us_histo, plastic_ds_histo, tuned_mc_sideband_signal, tuned_mc_sideband_Other, tuned_mc_sideband_regUS, tuned_mc_sideband_regDS, tmp_scale_factor, fixUS, fixDS, channelstr, onlyCalculateCV, writeSysCode ); //tune all universes as well

  //get directory of the scale factor and write scale factors histogram there 
  TDirectory *dir = scaleFactor->GetDirectory("");
  tmp_scale_factor->SetDirectory(dir);
 
  //tuned the scintillator background before adding it to the stacks 
  TObjArray hlistTuned;
  hlistTuned.Add(tuned_mc_sideband_signal);
  hlistTuned.Add(tuned_mc_sideband_Other);
  hlistTuned.Add(tuned_mc_sideband_regUS);
  hlistTuned.Add(tuned_mc_sideband_regDS);
 
/*        PlotUtils::HistFolio<PlotUtils::MnvH1D> hlistTuned(Form("tuned_grouped_%s_%s", targetString.c_str(), var.c_str()), tuned_mc_sideband_signal);
        hlistTuned.AddComponentHist("regDS", tuned_mc_sideband_regDS);
        hlistTuned.AddComponentHist("regUS", tuned_mc_sideband_regUS);
        hlistTuned.AddComponentHist("Other", tuned_mc_sideband_Other);
        hlistTuned.AddComponentHist("signal", tuned_mc_sideband_signal);
*/
      //PlotStacked(data_sideband, hlistTuned.GetHistArray(), dataMCScale, hlistTuned.GetName(), hlistTuned.GetName());

 
  //--- PLOT THE TUNED BACKGROUND --------
  // now let's draw the PLOT!
  //DrawMeSomeDataMCStackedHistos( mnvPlotter, data_sideband, hlistTuned, dataMCScale, canname, plotname, var, targetZ, makeDISPlots ); 
  
  //---- PLOT DATA/MC RATIO and CHI^2
  //PlotChi2Stat( mnvPlotter, data_sideband, hlistTuned, dataMCScale, var, targetZ, canname ); 
 
  TFile* tunedHists;
  if(makeDISPlots){ 
     tunedHists =  new TFile( Form("%s/TunedPlasticSidebandsDIS_%s_t%d_z%d_%s_%s.root", "./Plastics_TunedHist"/*outdir.c_str()*/, !fixUS? "US":"DS", targetID, targetZ, var.c_str(), playlist.c_str(), getenv("NUKECC_TAG") ), "recreate" );
   }
  else{
     tunedHists =  new TFile( Form("%s/TunedPlasticSidebands_%s_t%d_z%d_%s_%s.root", "./Plastics_TunedHist"/*outdir.c_str()*/, !fixUS? "US":"DS", targetID, targetZ, var.c_str(), playlist.c_str(), getenv("NUKECC_TAG") ), "recreate" );
   }

  tunedHists->cd();
  data_sideband->Write();
  tuned_mc_sideband_signal->Write();
  tuned_mc_sideband_Other->Write();
  tuned_mc_sideband_regUS->Write();
  tuned_mc_sideband_regDS->Write();

  tunedHists->Close();
}

MnvH1D*  TunePlasticBackground( MnvH1D *histo, string region, int targetZ, string playlist, string var ){
    
    TFile *scaleFactor = new TFile( Form("%s/CHScaleFactors_z%d_%s_%s.root", getenv("FILES"), targetZ, playlist.c_str(), getenv("NUKECC_TAG") ), "read" );
    assert(scaleFactor);
    
    MnvH1D* h_scaleFactors = (MnvH1D*)scaleFactor->Get( Form("scaleFactor_%s_%s", region.c_str(), var.c_str() ) ); //region should be "US" or "DS"
    cout << "name of scale factors = " << scaleFactor->GetName() << endl;
    int nonzerobin = h_scaleFactors->FindFirstBinAbove(0.);
    MnvH1D* tunedHisto = (MnvH1D*)histo->Clone("tunedHisto");
    MnvH1D* scaleHisto = (MnvH1D*)histo->Clone("scaleHisto");
    scaleHisto->Reset();
    cout<<" I cloned the histogram "<<endl;
    double cv_val = h_scaleFactors->GetBinContent(nonzerobin);
    cout<<"I've got the scale factor stuff for SB "<<region.c_str()<<"  the CV value is "<<cv_val<<endl;
    
    for(int i=0;i<scaleHisto->GetNbinsX()+2;i++) scaleHisto->SetBinContent(i,cv_val);
  
    tunedHisto->Multiply(tunedHisto,scaleHisto);
    
    return tunedHisto;
    
}

// ===========================================================
// you're entering background constraint zone 
// ============================================================

void  ConstrainScintBackground(MnvH1D* mnvh1d_data,
			       MnvH1D* mnvh1d_mc_signal,
			       MnvH1D* mnvh1d_mc_other,
			       MnvH1D* mnvh1d_mc_background_US,
			       MnvH1D* mnvh1d_mc_background_DS,
			       MnvH1D* mc_signal,
			       MnvH1D* mc_other,
			       MnvH1D* mc_background_US,
			       MnvH1D* mc_background_DS,
			       MnvH1D* tmp_scale_factor,
			       //double mc_scale,
			       bool fixUS,
			       bool fixDS,
			       string channel /*= "DIS"*/,
			       bool onlyCV /*=false*/,
			       bool writeScaleFactor /*=false*/ )
{
  
    //Get the CV histo
    TH1D* h1d_data_cv = new TH1D(mnvh1d_data->GetCVHistoWithStatError());
    TH1D* h1d_signal_cv = new TH1D(mnvh1d_mc_signal->GetCVHistoWithStatError());
    TH1D* h1d_other_cv = new TH1D(mnvh1d_mc_other->GetCVHistoWithStatError());
    TH1D* h1d_backgroundUS_cv = new TH1D(mnvh1d_mc_background_US->GetCVHistoWithStatError());
    TH1D* h1d_backgroundDS_cv = new TH1D(mnvh1d_mc_background_DS->GetCVHistoWithStatError());
    
    std::cout << "Total data: " << h1d_data_cv->Integral() << std::endl;
    
    vector<double> scales_cv, results_cv;
    // fit tells how much the MC should be scaled to match data
    results_cv = CalcScaleFactorMinimizer(h1d_data_cv,
                                          h1d_signal_cv,
                                          h1d_other_cv,
                                          h1d_backgroundUS_cv,
                                          h1d_backgroundDS_cv,
                                          //mc_scale,
                                          "CV",
                                          fixUS,
                                          fixDS);
    
    //scale only the cv here, universe histograms below
    //scale only the floated histogram;
    if( !fixUS && fixDS ){
        ((TH1D*)mc_background_US)->Scale(results_cv[0]);
        ((TH1D*)tmp_scale_factor)->Divide( ((TH1D*)tmp_scale_factor), ((TH1D*)tmp_scale_factor), results_cv[0], 1.0);
    }
    else if( !fixDS && fixUS ){
        ((TH1D*) mc_background_DS)->Scale(results_cv[0]);
        ((TH1D*)tmp_scale_factor)->Divide( ((TH1D*)tmp_scale_factor), ((TH1D*)tmp_scale_factor), results_cv[0], 1.0);
    }
    
    if( onlyCV ) return;

    
    // vertical error band
    std::vector<std::string> vertNames = mnvh1d_mc_signal->GetVertErrorBandNames();
cout<<__LINE__<<endl;
    for (std::vector<std::string>::iterator name = vertNames.begin();
         name != vertNames.end(); ++name) {
cout<<__LINE__<<endl;
        std::string errName = *name;
        const unsigned int n_universe = mnvh1d_mc_signal->GetVertErrorBand(errName)->GetNHists();
        
cout<<__LINE__<<endl;
        //Tune CV in each universe
        if( !fixUS && mnvh1d_mc_background_US->HasVertErrorBand(errName)) {
            
            MnvVertErrorBand *vertError = mc_background_US->GetVertErrorBand(errName);
            vertError->Scale(results_cv[0],"",false);
            MnvVertErrorBand *scaleVertError = tmp_scale_factor->GetVertErrorBand(errName);
            scaleVertError->TH1D::Divide( ((TH1D*)tmp_scale_factor), ((TH1D*)tmp_scale_factor), results_cv[0], 1.0);
            
        }
        if( !fixDS && mnvh1d_mc_background_DS->HasVertErrorBand(errName)) {
            
            MnvVertErrorBand *vertError = mc_background_DS->GetVertErrorBand(errName);
            vertError->Scale(results_cv[0],"",false);
            MnvVertErrorBand *scaleVertError = tmp_scale_factor->GetVertErrorBand(errName);
            scaleVertError->TH1D::Divide( ((TH1D*)tmp_scale_factor), ((TH1D*)tmp_scale_factor), results_cv[0], 1.0);
            
        }
        
        for (unsigned int i = 0; i < n_universe; ++i) {
            
            std::cout << "\t Uncertainty band name:   " << errName << " universe " << i << std::endl;
            std::cout << "\t\tTotal data:      " << h1d_data_cv->Integral() << std::endl;
            
            TH1D* h1d_mc_signal_universe        = GetVertErrorBandUniverseHist(mnvh1d_mc_signal,     errName, i);
            TH1D* h1d_mc_other_universe         = GetVertErrorBandUniverseHist(mnvh1d_mc_other,     errName, i);
            TH1D* h1d_mc_backgroundUS_universe  = GetVertErrorBandUniverseHist(mnvh1d_mc_background_US, errName, i);
            TH1D* h1d_mc_backgroundDS_universe  = GetVertErrorBandUniverseHist(mnvh1d_mc_background_DS, errName, i);
            
            vector<double> results_universe = CalcScaleFactorMinimizer( h1d_data_cv,
                                                                       h1d_mc_signal_universe,
                                                                       h1d_mc_other_universe,
                                                                       h1d_mc_backgroundUS_universe,
                                                                       h1d_mc_backgroundDS_universe,
                                                                       //mc_scale,
                                                                       errName,
                                                                       fixUS,
                                                                       fixDS  );
            
            assert(results_universe[0] > 0.0);
            
            // Check if this variable has systematic uncertainty with name 'errName'.
            // It can be false if the kinematic variable does not has the systematic uncertainty
            // Cannot just access the universe histogram and check if it is nullptr; it seems
            // that the universe histograms are not properly initialized to nullptr
            
            if( !fixUS && mnvh1d_mc_background_US->HasVertErrorBand(errName)) {
                
                TH1D* universe_histogram_US = GetVertErrorBandUniverseHist(mc_background_US, errName, i);
                TH1D* universe_scale_histo_US = GetVertErrorBandUniverseHist(tmp_scale_factor, errName, i);
                TH1D* universe_scalefactor_US = GetVertErrorBandUniverseHist(tmp_scale_factor, errName, i);
                
                //scale_histogram_bins(universe_histogram, scale_universe, boundary);
                universe_histogram_US->Scale(results_universe[0]);
                universe_scalefactor_US->Divide(universe_scale_histo_US,universe_scale_histo_US, results_universe[0], 1.0);
                
                std::cout.precision(8);
                std::cout << "\t\t " << errName << " " << i << " check scale US: " << results_universe[0] << std::endl;
                
            }// end of upstream
            
            if( !fixDS && mnvh1d_mc_background_DS->HasVertErrorBand(errName)) {
                
                TH1D* universe_histogram_DS = GetVertErrorBandUniverseHist(mc_background_DS, errName, i);
                TH1D* universe_scale_histo_DS = GetVertErrorBandUniverseHist(tmp_scale_factor, errName, i);
                TH1D* universe_scalefactor_DS = GetVertErrorBandUniverseHist(tmp_scale_factor, errName, i);
                
                //scale_histogram_bins(universe_histogram, scale_universe, boundary);
                universe_histogram_DS->Scale(results_universe[0]);
                universe_scalefactor_DS->Divide(universe_scale_histo_DS,universe_scale_histo_DS, results_universe[0], 1.0);
                
                std::cout.precision(8);
                std::cout << "\t\t  " << errName << " " << i << " check scale DS: " << results_universe[0] << std::endl;
            } //end of downstream
            
        }// end of universes
        
    }// end of vertical error band names
 
   
}

TH1D* GetVertErrorBandUniverseHist(MnvH1D* mnvh1d,
                                              const std::string& errName,
                                              int universe)
{
    //cout << "histo, errName, uni = " << mnvh1d << ", " << errName << ", " << universe << endl;
    return mnvh1d->GetVertErrorBand(errName)->GetHists().at(universe);
}


vector<double>  CalcScaleFactorMinimizer( TH1D* h1d_data,
					  TH1D* h1d_mc_signal,
					  TH1D* h1d_mc_other,
					  TH1D* h1d_mc_backgroundUS,
					  TH1D* h1d_mc_backgroundDS,
					  //double mc_scale, //already scaled
					  const std::string& errName,
					  bool fixUS /*=false*/,
					  bool fixDS /*=false*/ )
{
    
    m_histo_data = h1d_data;
    m_histo_sig_fix = h1d_mc_signal;
    m_histo_other_fix = h1d_mc_other;
    
    if( fixUS ){
        m_histo_plastic_fix = h1d_mc_backgroundUS;
        m_histo_plastic_float = h1d_mc_backgroundDS;
    }
    if( fixDS ){
        m_histo_plastic_fix = h1d_mc_backgroundDS;
        m_histo_plastic_float = h1d_mc_backgroundUS;
    }
    if( !fixDS && !fixUS ){
        m_histo_plastic_ds_float = h1d_mc_backgroundDS;
        m_histo_plastic_us_float = h1d_mc_backgroundUS;
    }
    
    // initialize stuff
    double scale_factor = 1.0;
    double minChi2 = -1.0;
    
    vector<double> results;
    results.push_back(scale_factor);
    results.push_back(minChi2);
    
    // Pick some numbers to be the tolerance and stuff
    ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2");
    fitter->SetMaxFunctionCalls(1000000);
    fitter->SetMaxIterations(100000);
    fitter->SetTolerance(0.01);
    fitter->SetPrintLevel(1);
    
    
    //BkgFitter::getChi2 has 1 parameter, minimizes chi2 for joint fit
    // initialize scale to 1.0 with step size of 0.001
    fitter->SetVariable( 0, "scale", 1.0, 0.001 );
    //ROOT::Math::Functor lf( this, & getChi2, 1 );
    ROOT::Math::Functor lf( & getChi2, 1 );
    ROOT::Math::Functor functor( lf, 1 );
    fitter->SetFunction( functor );
    
    // Go!
    fitter->Minimize();
    
    if( fitter->Status() != 0 ) {
        std::cout << "Womp womp womp... The background fit failed!" << std::endl;
        return results;
    }
    
    const double *bestfit = fitter->X();
    scale_factor = bestfit[0];
    minChi2 = getChi2( bestfit );
    results[0] = scale_factor;
    results[1] = minChi2;
    
    //cout << "Fitting using minimizer || scale factor, minChi2 = " << std::fixed << std::setprecision(3) << scale_factor << ", " << minChi2 << endl;
    return results;
    
}

double  getChi2( const double * par )
{
    
    double scale = par[0];
    
    // Sideband is bins from plane number 11 to 65
    // Assume no error on the MC, and use the stat errors on the data
    int binLow = m_histo_data->FindBin( 11. ); // same bins on both sets of histograms
    int binHigh = m_histo_data->FindBin( 65. );
    
    double chi2 = 0.0;
    // clone and scale the floating histograms
    TH1D * temp_plastic_float = new TH1D( *(TH1D*)m_histo_plastic_float->Clone("temp_plastic_float") );
    temp_plastic_float->Scale( scale );
    
    double tot_dat = 0.0;
    double tot_sig_fix = 0.0;
    double tot_other_fix = 0.0;
    double tot_plastic_fix = 0.0;
    double tot_plastic_float = 0.0;
    
    for( int b = binLow; b <= binHigh; ++b ) {
        double data = m_histo_data->GetBinContent(b);
        double err = m_histo_data->GetBinError(b);
        double sig_fix = m_histo_sig_fix->GetBinContent(b); // fixed part of MC, should be small
        double other_fix = m_histo_other_fix->GetBinContent(b); // fixed part of MC, should be small
        double plastic_fix = m_histo_plastic_fix->GetBinContent(b); // fixed part of MC, should be small
        double plastic_float = temp_plastic_float->GetBinContent(b); // let the background floats in the fitting
        
        if( data != 0.0 ){
            chi2 += pow( (plastic_float + plastic_fix + other_fix + sig_fix - data)/err, 2 ); // track bin chi2
            tot_dat += data;
            tot_sig_fix += sig_fix;
            tot_other_fix += other_fix;
            tot_plastic_fix += plastic_fix;
            tot_plastic_float += plastic_float;
        }
    }
    
    double unbinned_chi2 = pow( (tot_sig_fix+tot_other_fix+tot_plastic_fix+tot_plastic_float-tot_dat)/sqrt(tot_dat), 2 );
/*    
    bool floatUS = false;
    string histoName = m_histo_plastic_float->GetName();
    if( string::npos != histoName.find("plasticUS") ) floatUS = true;
    
    if( floatUS ){
        cout << "nominal event after fit ( data, signal, other, plastic US, plastic DS ) : " << std::fixed << std::setprecision(2) << tot_dat << " &  "
        << std::fixed << std::setprecision(2) <<  tot_sig_fix << " & "
        << std::fixed << std::setprecision(2) <<  tot_other_fix << " & "
        << std::fixed << std::setprecision(2) <<  tot_plastic_float << " & "
        << std::fixed << std::setprecision(2) <<  tot_plastic_fix << " & " << endl;
        cout << "nominal event after fit ( data, signal, other, plastic US, plastic DS ) : " << std::fixed << tot_dat << " & "
        << std::fixed << std::setprecision(2) <<  tot_sig_fix << " & "
        << std::fixed << std::setprecision(2) <<  tot_other_fix << " & "
        << std::fixed << std::setprecision(2) <<  tot_plastic_fix << " & "
        << std::fixed << std::setprecision(2) <<  tot_plastic_float << " & " << endl;

    }
*/    
    return unbinned_chi2;
    
}
/*
void AddHistosToAVectorUS( vector<MnvH1D*> &v_histos, vector<MnvH1D*> histos, int targetZ ){

  if( targetZ == 6 ) v_histos.push_back(planeDNN_US_C );
  if( targetZ == 26 ) v_histos.push_back(planeDNN_US_Fe );
  if( targetZ == 82 ) v_histos.push_back(planeDNN_US_Pb );

  v_histos.push_back(planeDNN_US_mat );
  v_histos.push_back(planeDNN_US_other );
  v_histos.push_back(planeDNN_US_regUS );
  v_histos.push_back(planeDNN_US_regDS );

}
void AddHistosToAVectorDS( vector<MnvH1D*> &v_histos, vector<MnvH1D*> histos, int targetZ ){

  if( targetZ == 6 ) v_histos.push_back(planeDNN_DS_C );
  if( targetZ == 26 ) v_histos.push_back(planeDNN_DS_Fe );
  if( targetZ == 82 ) v_histos.push_back(planeDNN_DS_Pb );

  v_histos.push_back(planeDNN_DS_mat );
  v_histos.push_back(planeDNN_DS_other );
  v_histos.push_back(planeDNN_DS_regUS );
  v_histos.push_back(planeDNN_DS_regDS );

}
*/
/*

void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetZ, TString canname ){
  
  //first declared some plotting parameters
  mnvPlotter.legend_offset_x = .03;
  mnvPlotter.width_xspace_per_letter = .33;
  
  string targetString;  
  if( targetZ == 6 )
    targetString = "Carbon";
  else if( targetZ == 26 )
    targetString = "Iron";
  if( targetZ == 82 )
    targetString = "Lead";
 
   
  //just get the CV for now 
  MnvH1D* data_cv = (MnvH1D*)dataHisto->Clone("mnvh1d_data_histo");
  cout << "data_cv = " << data_cv->Integral() << endl; 
  //get total MC
  MnvH1D* totalMC = (MnvH1D*)mchistos[0]->Clone("mnvh1d_mc_histo");
  totalMC->SetFillColor(0);
  totalMC->SetFillStyle(0);

  string cname = (string)canname; //"AfterTuning_..." plots
  if( targetZ == 6 && string::npos != cname.find("_US_") && string::npos != cname.find("_planeDNN_") ) totalMC->GetXaxis()->SetRangeUser(22., 28.);
  else if( targetZ == 6 && string::npos != cname.find("_DS_") && string::npos != cname.find("_planeDNN_") ) totalMC->GetXaxis()->SetRangeUser(34., 40.);
  
  for( int i = 1; i < mchistos.GetLast()+1; i++ ){
    MnvH1D* mc_histo_add = (MnvH1D*)mchistos[i]->Clone("mnvh1d_mc_histo_add");
    mc_histo_add->SetFillColor(0);
    mc_histo_add->SetFillStyle(0);
    totalMC->Add(mc_histo_add);
  }
  
//I do not know whats happening here...  
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
  cout<<"Label stat = "<<labelStat<<endl;
  
  string region = "Fiducial";
  if( string::npos != cname.find("_US_") ) region = "Upstream Plastic Background";
  if( string::npos != cname.find("_DS_") ) region = "Downstream Plastic Background";
  
  string status = "After Tuning";
  if( string::npos != cname.find("BeforeTuning_") ) status = "Before Tuning";

  if( writeChi2 && !tune_fiducial ) fprintf(Chi2Table, "%s %s %s %.2f \n", var.c_str(), region.c_str(), status.c_str(), chi2Sys/(Double_t)ndfStat ); 
  if( writeChi2 && tune_fiducial  ) fprintf(Chi2Table, "%s Fiducial %s, With Sys. %.2f \n", var.c_str(), status.c_str(), chi2Sys/(Double_t)ndfStat ); 

  double plotMin = -1., plotMax = -1.;
  GetMinAndMaxAxisRatio( data_cv, totalMC, dataMCScale, plotMin, plotMax );
  
  //plot data/mc with stat errors
  {
    TString cName = "BGPlotRatio_" + canname;
    TString Yaxis = "Data/MC (" + status + ")"; 
    TCanvas c( cName, cName, 1200, 800 );
    cout << "plotting " << cName << endl;
    mnvPlotter.DrawDataMCRatio( data_cv, totalMC, dataMCScale, true, true, plotMin, plotMax, Yaxis );
    mnvPlotter.AddPlotLabel( labelStat, .46, .875, .05 );
    
    if( WRITE_PRELIMINARY ) mnvPlotter.WritePreliminary("BC");
    
    if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("%s of %s", region.c_str(), targetString.c_str() ), .06);
    mnvPlotter.MultiPrint( &c );
  }
  
  {  
    TString cName = "FracError_" + canname;
    TCanvas can( cName, cName, 1200,800 );
    totalMC->SetFillColor(0);
    totalMC->SetFillStyle(0);
    totalMC->UseCurrentStyle();
    cout << "plotting " << cName << endl;
    totalMC->PopVertErrorBand(ErrorName::PlasticBG);
    //mnvPlotter.DrawDataMCErrorSummary( totalMC, data_cv, legLoc, false, true, 0.00001, false ); //use my colors and styles
    
    if(WRITE_PRELIMINARY)
      mnvPlotter.WritePreliminary(prelimLoc);
    
    if( ADD_HISTO_TITLES ){
      if( !tune_fiducial ) mnvPlotter.AddHistoTitle( Form("%s (%s), %s", region.c_str(), status.c_str(), targetString.c_str() ), 0.05 );
      else{
        if( string::npos != cname.find("_BGSubtracted_") ) mnvPlotter.AddHistoTitle( Form("Error Sources on Plastic Subtracted DIS Sample, %s", targetString.c_str() ), 0.04 );
        else mnvPlotter.AddHistoTitle( Form("Error Sources on Plastic Background (%s), %s", status.c_str(), targetString.c_str() ), 0.04 );
      }
    }
    mnvPlotter.MultiPrint( &can, can.GetName() );
  }
  
  //scale to data
  //totalMC->Scale(dataMCScale);
  {  
    TString cName = "FluxAllUniverse_" + canname;
    TCanvas can( cName, cName, 1200,800 );
    cout << "plotting " << cName << endl;
    //    mnvPlotter.DrawAllUniversesFlux( totalMC, data_cv, true ); //use my colors and styles
    data_cv->SetLineColor(kRed);
    data_cv->SetLineWidth(2);
    data_cv->Draw("hist same");
 
    if(WRITE_PRELIMINARY)
      mnvPlotter.WritePreliminary(prelimLoc);
    mnvPlotter.MultiPrint( &can, can.GetName() );
  }
   
}

TString GetPlotSuffix() const
{
    TString s = Form( "%s_%s%s", fPlaylist.c_str(), getenv("NUKECC_TAG"), GetVariationTag().Data() );
    return s;
}

*/
void DrawMeSomeDataMCStackedHistos( MnvPlotter mnvPlotter, MnvH1D* histo_data, TObjArray hlist, double scale, TString canname, string plotname, string var, int targetZ, const bool makeDISPlots ){
  
  //first declared some plotting parameters
  mnvPlotter.legend_offset_x = .03;
  mnvPlotter.width_xspace_per_letter = .33;
  
  //define max axis based on the sideband region 
  int opt = 0;
  if( string::npos != plotname.find(" Upstream ") ) opt = 1;
  if( string::npos != plotname.find(" Downstream ") ) opt = 2;
 
  cout << "canname ***** " << canname << endl; 
  TCanvas can( canname, canname, 1200,800 );
  //if( opt != 0 ) SetStackedHistosMaxAxis( mnvPlotter, targetZ, var, makeDISPlots, opt ); //set y axis of the plots for all the variables.
  mnvPlotter.DrawDataStackedMC( histo_data, &hlist, scale, legLoc, "Data", -1, -1, -1 ); //use my colors and styles
  
/*  mnvPlotter.WriteNorm("Abs-Normalized", infoX, infoY);
  mnvPlotter.WriteNorm("Stat. Errors Only", infoX, infoY-.04);
  if(WRITE_PRELIMINARY)
    mnvPlotter.WritePreliminary(prelimLoc);
  
  if( ADD_HISTO_TITLES )
    mnvPlotter.AddHistoTitle( plotname.c_str(), 0.042 );
*/  
  mnvPlotter.MultiPrint( &can, can.GetName() );

}


//Helper function to make the axis constant for specific kinematics and sepcific target
void SetStackedHistosMaxAxis( MnvPlotter &mnvPlotter, MnvH1D *histo, int targetZ, string var, const bool makeDISPlots, int opt ){
  
  //---------------------------------------------
  //  options:
  //---------------------------------------------
  //  1 = Upstream Sideband DIS
  //  2 = Downstream Sideband DIS
  //---------------------------------------------
  
  if( opt == 1 ){
    if( targetZ == 26 ) mnvPlotter.axis_maximum = 45000.0;                                                                                                                        
    if( targetZ == 82 ) mnvPlotter.axis_maximum = 50000.0;                                                                                                   
    if( targetZ == 6 ) mnvPlotter.axis_maximum = 25000.0;                                                                                                                       
    if( targetZ == 82 && var == "W" ) mnvPlotter.axis_maximum = 15000.0; 
    if( targetZ == 26 && var == "W" ) mnvPlotter.axis_maximum = 10000.0;
    if( targetZ == 6 && var == "W" ) mnvPlotter.axis_maximum = 5000.0; 
    if( targetZ == 82 && var == "ThetaMu" ) mnvPlotter.axis_maximum = 3000.0;
    if( targetZ == 26 && var == "ThetaMu" ) mnvPlotter.axis_maximum = 2000.0;
    if( targetZ == 6 && var == "ThetaMu" ) mnvPlotter.axis_maximum = 10000.0; 
    if( targetZ == 82 && var == "Ehad" ) mnvPlotter.axis_maximum = 2000.0;
    if( targetZ == 26 && var == "Ehad" ) mnvPlotter.axis_maximum = 15000.0;
    if( targetZ == 6 && var == "Ehad" ) mnvPlotter.axis_maximum = 10000.0;
    if( makeDISPlots ){
      if( targetZ == 82 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 3000.0;
      if( targetZ == 26 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 2000.0;
      if( targetZ == 6 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 10000.0;
    } else {
      if( targetZ == 82 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 30000.0;                                                                                                  
      if( targetZ == 26 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 25000.0;                                                                                                  
      if( targetZ == 6 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 15000.0;
    }
    if( targetZ == 82 && ( var == "Q2") ) mnvPlotter.axis_maximum = 1600.0;                                                                                                       
    if( targetZ == 26 && ( var == "Q2") ) mnvPlotter.axis_maximum = 1200.0;                                                                                                       
    if( targetZ == 6 && ( var == "Q2") ) mnvPlotter.axis_maximum = 800.0;                                                                                                         
    if( targetZ == 82 && (var == "x" || var == "y") ) mnvPlotter.axis_maximum = 800.0;                                                                                           
    if( targetZ == 26 && (var == "x" || var == "y") ) mnvPlotter.axis_maximum = 6000.0;                                                                                           
    if( targetZ == 6 && (var == "x" || var == "y") ) mnvPlotter.axis_maximum = 4000.0;                                                                                            
    if( targetZ == 82 && (var == "Emu" || var == "Emum") ) mnvPlotter.axis_maximum = 4000.0;                                                                                      
    if( targetZ == 26 && (var == "Emu" || var == "Emum") ) mnvPlotter.axis_maximum = 3000.0;                                                                                      
    if( targetZ == 6 && (var == "Emu" || var == "Emum") ) mnvPlotter.axis_maximum = 2000.0;
  }
  else if( opt == 2 ){
    if( targetZ == 26 ) mnvPlotter.axis_maximum = 45000.0; 
    if( targetZ == 82 ) mnvPlotter.axis_maximum = 50000.0; 
    if( targetZ == 6 ) mnvPlotter.axis_maximum = 25000.0; 
    if( targetZ == 26 && var == "W" ) mnvPlotter.axis_maximum = 15000.0; 
    if( targetZ == 82 && var == "W" ) mnvPlotter.axis_maximum = 10000.0; 
    if( targetZ == 6 && var == "W" ) mnvPlotter.axis_maximum = 5000.0; 
    if( targetZ == 26 && var == "ThetaMu" ) mnvPlotter.axis_maximum = 4000.0; 
    if( targetZ == 82 && var == "ThetaMu" ) mnvPlotter.axis_maximum = 250.0; 
    if( targetZ == 6 && var == "ThetaMu" ) mnvPlotter.axis_maximum = 15000.0; 
    if( targetZ == 26 && var == "Ehad" ) mnvPlotter.axis_maximum = 2000.0; 
    if( targetZ == 82 && var == "Ehad" ) mnvPlotter.axis_maximum = 15000.0; 
    if( targetZ == 6 && var == "Ehad" ) mnvPlotter.axis_maximum = 10000.0; 
    if( makeDISPlots ){ 
      if( targetZ == 82 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 3000.0; 
      if( targetZ == 26 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 2000.0; 
      if( targetZ == 6 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 10000.0;
    } else { 
      if( targetZ == 82 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 30000.0; 
      if( targetZ == 26 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 25000.0; 
      if( targetZ == 6 && ( var == "Enu" ) ) mnvPlotter.axis_maximum = 15000.0;
    }
    if( targetZ == 26 && ( var == "Q2") ) mnvPlotter.axis_maximum = 1600.0; 
    if( targetZ == 82 && ( var == "Q2") ) mnvPlotter.axis_maximum = 1200.0; 
    if( targetZ == 6 && ( var == "Q2") ) mnvPlotter.axis_maximum = 800.0; 
    if( targetZ == 26 && (var == "x" || var == "y") ) mnvPlotter.axis_maximum = 10000.0; 
    if( targetZ == 82 && (var == "x" || var == "y") ) mnvPlotter.axis_maximum = 6000.0; 
    if( targetZ == 6 && (var == "x" || var == "y") ) mnvPlotter.axis_maximum = 3000.0; 
    if( targetZ == 26 && (var == "Emu" || var == "Emum") ) mnvPlotter.axis_maximum = 4000.0; 
    if( targetZ == 82 && (var == "Emu" || var == "Emum") ) mnvPlotter.axis_maximum = 3000.0; 
    if( targetZ == 6 && (var == "Emu" || var == "Emum") ) mnvPlotter.axis_maximum = 2000.0;
    
  }
 
}
