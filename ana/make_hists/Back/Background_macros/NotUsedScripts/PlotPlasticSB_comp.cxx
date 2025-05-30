#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TArrayD.h" 
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
#include "include/CVUniverse.h"
#include "../../include/Variable_plasticSB.h"
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
#include "Cintex/Cintex.h"
#include "include/NukeCCUtilsNSF.h"

// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we
// shield it.
#ifndef __CINT__
#include "../../include/plotting_functions.h"
#endif
using namespace std;
using namespace NUKECC_ANA;
using namespace PlotUtils;

//----Helper Functions-----------------------------------------------
void PlotBGStuff(MnvPlotter mnvPlotter, vector<MnvH1D*> histos_mc, const string var, int targetID, int targetZ, double dataMCScale, bool plotUS, bool plotDS, string playlist, string Name);
//============================================================================================================================
// Main
int main(int argc, char * argv[]){
   ROOT::Cintex::Cintex::Enable();
   TH1::AddDirectory(false);
	
   if(argc==1){
     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
     std::cout<<"MACROS HELP:\n\n"<<
       "\t-./runEventLoop Path_to_Input_file Target_number Material_atomic_number Playlist\n\n"<<
       "\t-Path_to_Input_file\t =\t Path to the directory where the input ROOT file is called \n"<<
       "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" <<
       "\t-Playlist\t \t =\t eg. minervame1A"<< std::endl;
     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
     return 0;
   }

   string outdir=argv[1];
   int targetID = atoi(argv[2]);
   int targetZ = atoi(argv[3]);
   const string playlist= argv[4];
   
   ROOT::Cintex::Cintex::Enable();
   TH1::AddDirectory(false);
   
   bool RunCodeWithSystematics = false;
   TString histFileName, SF_file;

   cout<<"******************************************************************************************"<<endl;
   cout<<"                   I am making plots for both before and after tuning!!                   "<<endl;
   cout<<"******************************************************************************************"<<endl;

   if(RunCodeWithSystematics){
       histFileName = Form("%s/Hists_PlasticBackgd_with_SYS_t%d_z%02d_Nu_%s.root", outdir.c_str(), targetID, targetZ, getenv("NUKECC_TAG") ); 
     }
   else{
       histFileName = Form("%s/Hists_PlasticBackgd_without_SYS_t%d_z%02d_Nu_%s.root", outdir.c_str(), targetID, targetZ, getenv("NUKECC_TAG"));
     } 
   
   cout<<histFileName<<endl;

   MnvPlotter mnvPlotter(kNukeCCStyle);
   //MnvPlotter mnvPlotter(kCCQENuStyle);

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
   vars.push_back("Emu");

   double dataMCScale = 0.433;
   TFile *f1 = new TFile( histFileName,"read" );

   vector< MnvH1D*> hists_US_data, hists_US_mat, hists_US_other, hists_US_regUS, hists_US_regDS;
   vector< MnvH1D*> hists_DS_data, hists_DS_mat, hists_DS_other, hists_DS_regUS, hists_DS_regDS;
   vector< MnvH1D*> histos_mc_US, histos_mc_DS;
   vector< MnvH1D*> scale_factor_US, scale_factor_DS;
   vector< MnvH1D*> sf_US, sf_DS;

for( vector<string>::iterator i = vars.begin(); i != vars.end(); ++i ){
   const string& var = *i;


   f1->cd();
   hists_US_data.push_back(    (MnvH1D*)f1->Get(Form("selected_data_reco_US_%s", var.c_str())));
   hists_US_mat.push_back(     (MnvH1D*)f1->Get(Form("US_t2_z26_%s_Fe", var.c_str())));
   hists_US_other.push_back(   (MnvH1D*)f1->Get(Form("US_t2_z26_%s_Other", var.c_str())));
   hists_US_regUS.push_back(   (MnvH1D*)f1->Get(Form("US_t2_z26_%s_regUS", var.c_str())));
   hists_US_regDS.push_back(   (MnvH1D*)f1->Get(Form("US_t2_z26_%s_regDS", var.c_str())));
  
   hists_DS_data.push_back(  (MnvH1D*)f1->Get(Form("selected_data_reco_DS_%s", var.c_str())));
   hists_DS_mat.push_back(     (MnvH1D*)f1->Get(Form("DS_t2_z26_%s_Fe", var.c_str())));
   hists_DS_other.push_back(   (MnvH1D*)f1->Get(Form("DS_t2_z26_%s_Other", var.c_str())));
   hists_DS_regUS.push_back(   (MnvH1D*)f1->Get(Form("DS_t2_z26_%s_regUS", var.c_str())));
   hists_DS_regDS.push_back(   (MnvH1D*)f1->Get(Form("DS_t2_z26_%s_regDS", var.c_str())));
    }

for( unsigned int i = 0; i != vars.size(); ++i ){
  const string& var = vars[i];
  cout<<" variable in the loop = "<<var<<endl;

  histos_mc_US.push_back( hists_US_data[i] );
  histos_mc_US.push_back( hists_US_mat[i] );
  histos_mc_US.push_back( hists_US_other[i] );
  histos_mc_US.push_back( hists_US_regUS[i] );
  histos_mc_US.push_back( hists_US_regDS[i] );
      
  histos_mc_DS.push_back( hists_DS_data[i] );
  histos_mc_DS.push_back( hists_DS_mat[i] );
  histos_mc_DS.push_back( hists_DS_other[i] );
  histos_mc_DS.push_back( hists_DS_regUS[i] );
  histos_mc_DS.push_back( hists_DS_regDS[i] );
 
  PlotBGStuff(mnvPlotter, histos_mc_US, vars[i], targetID, targetZ, dataMCScale, true, false, playlist, "Untuned" );  //Filling US sidebands
  PlotBGStuff(mnvPlotter, histos_mc_DS, vars[i], targetID, targetZ, dataMCScale, false, true, playlist, "UnTuned");   //Filling DS sidebands

    //clear vector before going to the next variable 
    histos_mc_US.clear();
    histos_mc_DS.clear();
 
  }

// Now I am going to make tuned plots!!!! 
  vector< MnvH1D*> hists_US_mat_tuned, hists_US_other_tuned, hists_US_regUS_tuned, hists_US_regDS_tuned;
  vector< MnvH1D*> hists_DS_mat_tuned, hists_DS_other_tuned, hists_DS_regUS_tuned, hists_DS_regDS_tuned;
  vector< MnvH1D*> histos_mc_US_tuned, histos_mc_DS_tuned;
  vector< MnvH1D*> hists_US_data_tuned, hists_DS_data_tuned;
    bool plotUS = true;
    bool plotDS = true;
  for( vector<string>::iterator i = vars.begin(); i != vars.end(); ++i ){  
    const string& var = *i;
    TString histFileNameUS = Form("%s/TunedPlasticSidebands_US_t%d_z%02d_%s_%s.root", outdir.c_str(), targetID, targetZ, var.c_str(), playlist.c_str());
    TString histFileNameDS = Form("%s/TunedPlasticSidebands_DS_t%d_z%02d_%s_%s.root", outdir.c_str(), targetID, targetZ, var.c_str(), playlist.c_str());
 
    TFile *f_us = new TFile( histFileNameUS,"read" );
    TFile *f_ds = new TFile( histFileNameDS,"read" );

    f_us->cd();
    hists_US_data_tuned.push_back(  (MnvH1D*)f_us->Get("tuned_data"));
    hists_US_mat_tuned.push_back(   (MnvH1D*)f_us->Get("tuned_mc_signal"));
    hists_US_other_tuned.push_back( (MnvH1D*)f_us->Get("tuned_mc_other"));
    hists_US_regUS_tuned.push_back( (MnvH1D*)f_us->Get("tuned_mc_us"));
    hists_US_regDS_tuned.push_back( (MnvH1D*)f_us->Get("tuned_mc_ds"));
   
    f_ds->cd();
    hists_DS_data_tuned.push_back(  (MnvH1D*)f_ds->Get("tuned_data"));
    hists_DS_mat_tuned.push_back(   (MnvH1D*)f_ds->Get("tuned_mc_signal"));
    hists_DS_other_tuned.push_back( (MnvH1D*)f_ds->Get("tuned_mc_other"));
    hists_DS_regUS_tuned.push_back( (MnvH1D*)f_ds->Get("tuned_mc_us"));
    hists_DS_regDS_tuned.push_back( (MnvH1D*)f_ds->Get("tuned_mc_ds"));
 
  }


  for( unsigned int i = 0; i != vars.size(); ++i ){  
    const string& var = vars[i];
    cout<<" variable in the after tuning loop = "<<var<<endl;
  
    histos_mc_US_tuned.push_back( hists_US_data_tuned[i] );
    histos_mc_US_tuned.push_back( hists_US_mat_tuned[i] );
    histos_mc_US_tuned.push_back( hists_US_other_tuned[i] );
    histos_mc_US_tuned.push_back( hists_US_regUS_tuned[i] );
    histos_mc_US_tuned.push_back( hists_US_regDS_tuned[i] );
      
    histos_mc_DS_tuned.push_back( hists_DS_data_tuned[i] );
    histos_mc_DS_tuned.push_back( hists_DS_mat_tuned[i] );
    histos_mc_DS_tuned.push_back( hists_DS_other_tuned[i] );
    histos_mc_DS_tuned.push_back( hists_DS_regUS_tuned[i] );
    histos_mc_DS_tuned.push_back( hists_DS_regDS_tuned[i] );

  
    PlotBGStuff(mnvPlotter, histos_mc_US_tuned, vars[i], targetID, targetZ, dataMCScale, true, false, playlist, "Tuned" );   //Filling US sidebands
    PlotBGStuff(mnvPlotter, histos_mc_DS_tuned, vars[i], targetID, targetZ, dataMCScale, false, true, playlist, "Tuned" );   //Filling DS sidebands
    
    //clear vector before going to the next variable 
      histos_mc_US_tuned.clear();
      histos_mc_DS_tuned.clear();
            
    }
}// End of Main

void PlotBGStuff(MnvPlotter mnvPlotter, vector<MnvH1D*> histos_mc, const string var, int targetID, int targetZ, double dataMCScale, bool plotUS, bool plotDS, string playlist, string Name){
 
  MnvH1D* data_sideband       = (MnvH1D*)histos_mc[0]->Clone("data");
  MnvH1D* mc_sideband_signal  = (MnvH1D*)histos_mc[1]->Clone("mc_signal");
  MnvH1D* mc_sideband_Other   = (MnvH1D*)histos_mc[2]->Clone("mc_other");
  MnvH1D* mc_sideband_regUS   = (MnvH1D*)histos_mc[3]->Clone("mc_us");
  MnvH1D* mc_sideband_regDS   = (MnvH1D*)histos_mc[4]->Clone("mc_us");

  TObjArray hlist;
  hlist.Add(mc_sideband_signal);
  hlist.Add(mc_sideband_Other);
  hlist.Add(mc_sideband_regUS);
  hlist.Add(mc_sideband_regDS);


  string name;
  if(plotUS) name = Form("US_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
  if(plotDS) name = Form("DS_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
 

  // Plotting sidebands  
  PlotStacked(data_sideband, hlist, dataMCScale, var.c_str(), name, "");


}

