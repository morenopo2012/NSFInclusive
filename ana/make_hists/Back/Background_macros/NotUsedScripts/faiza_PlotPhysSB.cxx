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
#include "../../include/Variable_physicsSB.h"
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
void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name);
void GetMinAndMaxAxisRatio( TH1D* histData, TH1D* histMC, double dataMCScale, double &plotMin, double &plotMax );

void PlotBGRatioStuff(MnvPlotter mnvPlotter, MnvH1D* histos_mc, MnvH1D* histos_mc_tuned, const string var, int targetID, int targetZ, double dataMCScale, double dataPOT, double mcPOT, bool plotUS, bool plotDS, string playlist, string Name);
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
   
   //bool RunCodeWithSystematics = true;
   bool RunCodeWithSystematics = false;
   TString histFileName;

   cout<<"******************************************************************************************"<<endl;
   cout<<"                   I am making plots for both before and after tuning!!                   "<<endl;
   cout<<"******************************************************************************************"<<endl;

   if(RunCodeWithSystematics){
       histFileName = Form("%s/Hists_PhysicsBackgd_with_SYS_FullDet_q2WfromBranch_ME1A_targetsCombined_t%d_z%02d_Nu_%s.root", outdir.c_str(), targetID, targetZ, getenv("NUKECC_TAG") ); 
     }
   else{
       histFileName = Form("%s/Hists_PhysicsBackgd_without_SYS_t%d_z%02d_Nu_%s.root", outdir.c_str(), targetID, targetZ, getenv("NUKECC_TAG"));
     } 
   
   cout<<histFileName<<endl;

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
   vars.push_back("Emu");
   vars.push_back("Enu");
   vars.push_back("Ehad");

   TFile *f1 = new TFile( histFileName,"read" );

   TParameter<double> *mcPOT = (TParameter<double>*)f1->Get("MCPOT");
   TParameter<double> *dataPOT = (TParameter<double>*)f1->Get("DataPOT");
   double mcpot = mcPOT->GetVal();
   double datapot = dataPOT->GetVal(); 
   double dataMCScale = datapot/mcpot;
   cout<<"MCPOT = "<<mcpot<<"DataPOT = "<< datapot << "Scale = " << dataMCScale<<endl;
   
   TObjArray hlistUntunedTrans, hlistUntunedContin, hlisttunedTrans, hlisttunedContin;


   vector< MnvH1D*> hists_Trans_data, hists_Trans_in_Trans, hists_Contin_in_Trans, hists_Signal_in_Trans;
   vector< MnvH1D*> hists_Contin_data, hists_Trans_in_Contin, hists_Contin_in_Contin, hists_Signal_in_Contin;
   vector< MnvH1D*> histos_mc_trans, histos_mc_contin;

for( vector<string>::iterator i = vars.begin(); i != vars.end(); ++i ){
   const string& var = *i;

   f1->cd();
   hists_Trans_data.push_back(      (MnvH1D*)f1->Get(Form("selected_data_reco_trans_%s", var.c_str())));
   hists_Trans_in_Trans.push_back(  (MnvH1D*)f1->Get(Form("hists_trans_in_trans_%s", var.c_str())));
   hists_Contin_in_Trans.push_back( (MnvH1D*)f1->Get(Form("hists_contin_in_trans_%s", var.c_str())));
   hists_Signal_in_Trans.push_back( (MnvH1D*)f1->Get(Form("hists_signal_in_trans_%s",  var.c_str())));
  
   hists_Contin_data.push_back(      (MnvH1D*)f1->Get(Form("selected_data_reco_contin_%s", var.c_str())));
   hists_Trans_in_Contin.push_back(  (MnvH1D*)f1->Get(Form("hists_trans_in_contin_%s", var.c_str())));
   hists_Contin_in_Contin.push_back( (MnvH1D*)f1->Get(Form("hists_contin_in_contin_%s", var.c_str())));
   hists_Signal_in_Contin.push_back( (MnvH1D*)f1->Get(Form("hists_signal_in_contin_%s",  var.c_str())));
 
   MnvH1D *hists_Trans_data      = (MnvH1D*)f1->Get(Form("selected_data_reco_trans_%s", var.c_str()));
   MnvH1D *hists_Trans_in_Trans  = (MnvH1D*)f1->Get(Form("hists_trans_in_trans_%s", var.c_str()));
   MnvH1D *hists_Contin_in_Trans = (MnvH1D*)f1->Get(Form("hists_contin_in_trans_%s", var.c_str()));
   MnvH1D *hists_Signal_in_Trans = (MnvH1D*)f1->Get(Form("hists_signal_in_trans_%s",  var.c_str()));

   //save them into TObjArray
   TObjArray hlistUntunedTrans;
   hlistUntunedTrans.Add(hists_Trans_in_Trans);
   hlistUntunedTrans.Add(hists_Contin_in_Trans);
   hlistUntunedTrans.Add(hists_Signal_in_Trans);

   MnvH1D *hists_Contin_data      = (MnvH1D*)f1->Get(Form("selected_data_reco_contin_%s", var.c_str()));
   MnvH1D *hists_Trans_in_Contin  = (MnvH1D*)f1->Get(Form("hists_trans_in_contin_%s", var.c_str()));
   MnvH1D *hists_Contin_in_Contin = (MnvH1D*)f1->Get(Form("hists_contin_in_contin_%s", var.c_str()));
   MnvH1D *hists_Signal_in_Contin = (MnvH1D*)f1->Get(Form("hists_signal_in_contin_%s",  var.c_str()));  
  
   //save them into TObjArray
   TObjArray hlistUntunedContin;
   hlistUntunedContin.Add(hists_Trans_in_Contin);
   hlistUntunedContin.Add(hists_Contin_in_Contin);
   hlistUntunedContin.Add(hists_Signal_in_Contin);

  PlotChi2Stat( mnvPlotter, hists_Trans_data, hlistUntunedTrans, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "Trans_Untuned");
  PlotChi2Stat( mnvPlotter, hists_Contin_data, hlistUntunedContin, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "Contin_Untuned");

   }

for( unsigned int i = 0; i != vars.size(); ++i ){
  const string& var = vars[i];
  cout<<" variable in the loop = "<<var<<endl;

  histos_mc_trans.push_back( hists_Trans_data[i] );
  histos_mc_trans.push_back( hists_Signal_in_Trans[i] );
  histos_mc_trans.push_back( hists_Trans_in_Trans[i] );
  histos_mc_trans.push_back( hists_Contin_in_Trans[i] );
      
  histos_mc_contin.push_back( hists_Contin_data[i] );
  histos_mc_contin.push_back( hists_Signal_in_Contin[i] );
  histos_mc_contin.push_back( hists_Trans_in_Contin[i] );
  histos_mc_contin.push_back( hists_Contin_in_Contin[i] );
 

  PlotBGStuff(mnvPlotter, histos_mc_trans, vars[i], targetID, targetZ, dataMCScale, true, false, playlist, "Untuned" );  //Filling trans sidebands
  PlotBGStuff(mnvPlotter, histos_mc_contin, vars[i], targetID, targetZ, dataMCScale, false, true, playlist, "UnTuned");   //Filling contin sidebands

    //clear vector before going to the next variable 
    histos_mc_trans.clear();
    histos_mc_contin.clear();
 
  }

// Now I am going to make tuned plots!!!! 
   vector< MnvH1D*> hists_Signal_in_Trans_tuned, hists_Trans_in_Trans_tuned, hists_Contin_in_Trans_tuned;
   vector< MnvH1D*> hists_Signal_in_Contin_tuned, hists_Trans_in_Contin_tuned, hists_Contin_in_Contin_tuned;
   vector< MnvH1D*> histos_mc_trans_tuned, histos_mc_contin_tuned;
   vector< MnvH1D*> hists_Trans_data_tuned, hists_Contin_data_tuned;

  for( vector<string>::iterator i = vars.begin(); i != vars.end(); ++i ){  
    const string& var = *i;

    TString histFilename = Form("%s/TunedPhysicsSidebands_t%d_z%02d_%s_%s_%s.root", outdir.c_str(), targetID, targetZ, var.c_str(), playlist.c_str(), getenv("NUKECC_TAG"));
    TFile *f_tuned = new TFile( histFilename,"read" );
    cout<<"After Tuning filename = "<<histFilename<<endl;

    f_tuned->cd();
    hists_Trans_data_tuned.push_back(       (MnvH1D*)f_tuned->Get(Form("selected_data_reco_trans_%s", var.c_str())));
    hists_Signal_in_Trans_tuned.push_back(  (MnvH1D*)f_tuned->Get("sig_in_trans"));
    hists_Trans_in_Trans_tuned.push_back(   (MnvH1D*)f_tuned->Get("trans_in_trans"));
    hists_Contin_in_Trans_tuned.push_back(  (MnvH1D*)f_tuned->Get("contin_in_trans"));
 
    hists_Contin_data_tuned.push_back(       (MnvH1D*)f_tuned->Get(Form("selected_data_reco_contin_%s", var.c_str())));
    hists_Signal_in_Contin_tuned.push_back(  (MnvH1D*)f_tuned->Get("sig_in_contin"));
    hists_Trans_in_Contin_tuned.push_back(   (MnvH1D*)f_tuned->Get("trans_in_contin"));
    hists_Contin_in_Contin_tuned.push_back(  (MnvH1D*)f_tuned->Get("contin_in_contin"));

    f_tuned->cd();
    MnvH1D *hists_Trans_data_tuned      = (MnvH1D*)f_tuned->Get(Form("selected_data_reco_trans_%s", var.c_str()));
    MnvH1D *hists_Signal_in_Trans_tuned = (MnvH1D*)f_tuned->Get("sig_in_trans");
    MnvH1D *hists_Trans_in_Trans_tuned  = (MnvH1D*)f_tuned->Get("trans_in_trans");
    MnvH1D *hists_Contin_in_Trans_tuned = (MnvH1D*)f_tuned->Get("contin_in_trans");
  
    //save them into TObjArray
    TObjArray hlisttunedTrans;
    hlisttunedTrans.Add(hists_Trans_in_Trans_tuned);
    hlisttunedTrans.Add(hists_Contin_in_Trans_tuned);
    hlisttunedTrans.Add(hists_Signal_in_Trans_tuned);
 
    MnvH1D *hists_Contin_data_tuned      = (MnvH1D*)f_tuned->Get(Form("selected_data_reco_contin_%s", var.c_str()));
    MnvH1D *hists_Signal_in_Contin_tuned = (MnvH1D*)f_tuned->Get("sig_in_contin");
    MnvH1D *hists_Trans_in_Contin_tuned  = (MnvH1D*)f_tuned->Get("trans_in_contin");
    MnvH1D *hists_Contin_in_Contin_tuned = (MnvH1D*)f_tuned->Get("contin_in_contin");

    //save them into TObjArray
    TObjArray hlisttunedContin;
    hlisttunedContin.Add(hists_Trans_in_Contin_tuned);
    hlisttunedContin.Add(hists_Contin_in_Contin_tuned);
    hlisttunedContin.Add(hists_Signal_in_Contin_tuned);

    PlotChi2Stat( mnvPlotter, hists_Trans_data_tuned, hlisttunedTrans, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "Trans_tuned");

    PlotChi2Stat( mnvPlotter, hists_Contin_data_tuned, hlisttunedContin, dataMCScale, var.c_str(), targetID, targetZ, true, false, playlist, "Contin_tuned");
 
}
  
  for( unsigned int i = 0; i != vars.size(); ++i ){  
    const string& var = vars[i];
    cout<<" variable in the after tuning loop = "<<var<<endl;
 
 
    histos_mc_trans_tuned.push_back( hists_Trans_data_tuned[i] );
    histos_mc_trans_tuned.push_back( hists_Signal_in_Trans_tuned[i] );
    histos_mc_trans_tuned.push_back( hists_Trans_in_Trans_tuned[i] );
    histos_mc_trans_tuned.push_back( hists_Contin_in_Trans_tuned[i] );
      
    histos_mc_contin_tuned.push_back( hists_Contin_data_tuned[i] );
    histos_mc_contin_tuned.push_back( hists_Signal_in_Contin_tuned[i] );
    histos_mc_contin_tuned.push_back( hists_Trans_in_Contin_tuned[i] );
    histos_mc_contin_tuned.push_back( hists_Contin_in_Contin_tuned[i] );
  
    PlotBGStuff(mnvPlotter, histos_mc_trans_tuned, vars[i], targetID, targetZ, 1.0, true, false, playlist, "Tuned" );   //Filling US sidebands
    PlotBGStuff(mnvPlotter, histos_mc_contin_tuned, vars[i], targetID, targetZ,1.0, false, true, playlist, "Tuned" );   //Filling DS sidebands
    
    //clear vector before going to the next variable 
      histos_mc_trans_tuned.clear();
      histos_mc_contin_tuned.clear();
            
    }

}// End of Main

void PlotBGStuff(MnvPlotter mnvPlotter, vector<MnvH1D*> histos_mc, const string var, int targetID, int targetZ, double dataMCScale, bool plotUS, bool plotDS, string playlist, string Name){
 
  MnvH1D* data_sideband       = (MnvH1D*)histos_mc[0]->Clone("data");
  MnvH1D* mc_sideband_signal  = (MnvH1D*)histos_mc[1]->Clone("mc_signal");
  MnvH1D* mc_sideband_trans   = (MnvH1D*)histos_mc[2]->Clone("mc_trans");
  MnvH1D* mc_sideband_contin  = (MnvH1D*)histos_mc[3]->Clone("mc_contin");

  string suffix = "Faiza_is_trying_again";
  PlotUtils::HistFolio<PlotUtils::MnvH1D> hlist_grouped(Form("%s_%s", suffix.c_str(), var.c_str()), mc_sideband_signal);
    hlist_grouped.AddComponentHist("Signal", mc_sideband_signal);
    hlist_grouped.AddComponentHist("Trans",  mc_sideband_trans);
    hlist_grouped.AddComponentHist("Contin", mc_sideband_contin);

  string name = "Faiza_is_trying";
  if(plotUS) name = Form("Trans_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
  if(plotDS) name = Form("Contin_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
 

  // Plotting sidebands  
  PlotStacked(data_sideband, hlist_grouped.GetHistArray(), dataMCScale, var.c_str(), name, "");


}

void PlotBGRatioStuff(MnvPlotter mnvPlotter, MnvH1D* histos_mc, MnvH1D* histos_mc_tuned, const string var, int targetID, int targetZ, double dataMCScale, double dataPOT, double mcPOT, bool plotUS, bool plotDS, string playlist, string Name){

  string name;
  if(plotUS) name = Form("Ratio_Trans_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());
  if(plotDS) name = Form("Ratio_Contin_%s_%s_t%d_z%02d_%s", Name.c_str(), var.c_str(), targetID, targetZ, playlist.c_str());

  PlotRatio(histos_mc, histos_mc_tuned, dataMCScale, dataPOT, mcPOT, var.c_str(), name, "");
  //PlotRatio(histos_mc_tuned, histos_mc, dataMCScale, var.c_str(), name, "");


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
    //mnvPlotter.DrawDataMCRatio( data_cv, totalMC, dataMCScale, true, true, plotMin, plotMax, Yaxis );
    mnvPlotter.DrawDataMCRatio( data_cv, totalMC, dataMCScale, true, true, 0., 5., x_label.c_str(), Yaxis );
    mnvPlotter.AddPlotLabel( labelStat, .46, .875, .05 );
    
    if( WRITE_PRELIMINARY ) mnvPlotter.WritePreliminary("BC");
    
    //if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("%s of %s", region.c_str(), targetString.c_str() ), .06);
    mnvPlotter.MultiPrint( &c );
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




