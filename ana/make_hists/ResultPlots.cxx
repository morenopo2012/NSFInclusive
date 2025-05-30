#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TArrayD.h" 
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "../include/Variable_plasticSB.h"
#include "../include/Variable_physicsSB.h"
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
#include "../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"

// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we
// shield it.
#ifndef __CINT__
#include "../include/plotting_functions.h"
#endif
using namespace std;
using namespace NUKECC_ANA;
using namespace PlotUtils;

//----Helper Functions-----------------------------------------------
void PlotBGStuff(MnvPlotter mnvPlotter, vector<MnvH1D*> histos_mc, const string var, int targetID, int targetZ, double dataMCScale, bool plotUS, bool plotDS, string playlist, string Name);
void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name);
void GetMinAndMaxAxisRatio( TH1D* histData, TH1D* histMC, double dataMCScale, double &plotMin, double &plotMax );

MnvH1D *SumMaterialHists( const std::string fName, const std::string hName, const std::string var, FileType::t_FileType type, int targetZ, string playlist, bool subtractNonDIS /*=false*/ );
void PlotBGRatioStuff(MnvPlotter mnvPlotter, MnvH1D* histos_mc, MnvH1D* histos_mc_tuned, const string var, int targetID, int targetZ, double dataMCScale, double dataPOT, double mcPOT, bool plotUS, bool plotDS, string playlist, string Name);

//============================================================================================================================
// Main
int main(int argc, char * argv[]){
   //ROOT::Cintex::Cintex::Enable();
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
   
   //ROOT::Cintex::Cintex::Enable();
   TH1::AddDirectory(false);
   
   const bool sumMaterial = (targetID == 3 || targetID == 14);
   bool RunCodeWithSystematics = true;
   //bool RunCodeWithSystematics = false;
   TString histFileName, effFileName, eventfile;

   cout<<"******************************************************************************************"<<endl;
   cout<<"                   I am making plots for both before and after tuning!!                   "<<endl;
   cout<<"******************************************************************************************"<<endl;

   if(RunCodeWithSystematics){
       histFileName = Form("%s/Unfolded_t%d_z%02d_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str() ); 
       effFileName = Form("%s/Hists_Efficiency_t%d_z%02d_Nu.root", outdir.c_str(), targetID, targetZ ); 
       eventfile = Form("%s/Hists_EventSelection_t%d_z%02d_Nu.root", outdir.c_str(), targetID, targetZ ); 
     }
   else{
       //histFileName = Form("%s/Hists_PhysicsBackgd_without_SYS_t%d_z%02d_Nu_%s.root", outdir.c_str(), targetID, targetZ, getenv("NUKECC_TAG"));
       histFileName = Form("%s/Hists_PhysicsBackgd_without_SYS_t%d_z%02d_Nu.root", outdir.c_str(), targetID, targetZ);
     } 
   
   cout<<histFileName<<endl;

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

   vector<string> vars;
   vars.clear();
//   vars.push_back("Emu");
   vars.push_back("Enu");
//   vars.push_back("Ehad");
   vars.push_back("x");
//   vars.push_back("y");
  

   
 
    std::ofstream mcxBjTable;
    
    if(targetID < 10)
        mcxBjTable.open(Form("Event_table_z%02d.tex", targetZ), std::ofstream::out);
    
    else
        mcxBjTable.open(Form("Event_table_z%02d.tex", -1), std::ofstream::out); 
   
  string targetString;  
  if( targetZ == 6 )
    targetString = "Carbon";
  else if( targetZ == 26 )
    targetString = "Iron of target 1";
  if( targetZ == 82 )
    targetString = "Lead";

   //const string targetString = targetID < 10 ?

   const string suffix = Form( "t%02d_z%02d", targetID, targetZ );
   
   HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(playlist);

   TFile *f0 = new TFile( eventfile,"read" );
   assert(f0);
   mnvPlotter.ApplyAxisStyle( f0 );
   
   TFile *f1 = new TFile( histFileName,"read" );
   assert(f1);
   mnvPlotter.ApplyAxisStyle( f1 );
 
   TFile *Eff = new TFile(effFileName);
   assert(Eff);
   mnvPlotter.ApplyAxisStyle(Eff);
 
  const double nnuc_MC = TargetUtils::Get().GetPassiveTargetNNucleons( targetID, targetZ, false );
  const double nnuc_D = TargetUtils::Get().GetPassiveTargetNNucleons( targetID, targetZ, false );
     
   TParameter<double> *mcPOT = (TParameter<double>*)f0->Get("MCPOT");
   TParameter<double> *dataPOT = (TParameter<double>*)f0->Get("DataPOT");
   double mcpot = mcPOT->GetVal();
   double datapot = dataPOT->GetVal(); 
   double dataMCScale = datapot/mcpot;
   cout<<"MCPOT = "<<mcpot<<"DataPOT = "<< datapot << "Scale = " << dataMCScale<<endl;
   
   TObjArray hlistUntunedTrans, hlistUntunedContin, hlisttunedTrans, hlisttunedContin;


   std::vector<MnvH1D*> dHists, dHists_unfolded;
   std::vector<MnvH1D*> mcHists, mcHists_unfolded;
   vector< MnvH1D*> hists_Trans_data, hists_Trans_in_Trans, hists_Contin_in_Trans, hists_Signal_in_Trans;
   vector< MnvH1D*> hists_Contin_data, hists_Trans_in_Contin, hists_Contin_in_Contin, hists_Signal_in_Contin;
   vector< MnvH1D*> histos_mc_trans, histos_mc_contin, data_sideband, data_sideband_Contin;
   std::vector<MnvH1D*> effHists;

  //Get the cross section hists
  std::vector<MnvH1D*> xsecHists;
  std::vector<MnvH1D*> evtrateHists;
  //TFile *fRawXSec = new TFile( utils->GetXSecFileName().Data() );
  //TFile *fRawXSec = new TFile( Form("%s/CrossSection_z%d_%s.root", outdir.c_str(), targetZ, playlist.c_str() ), "recreate" );
    cout<< "I'm looking for this GENIE file " << endl;

  
  
  TFile *fRawXSec = new TFile( utils->GetXSecFileName().Data() );
    cout<< "I'm looking for this GENIE file "<< utils->GetXSecFileName().Data() << endl;
assert( fRawXSec );
   const bool useBG = targetID < 10;
  
   for( vector<string>::iterator ivar = vars.begin(); ivar != vars.end(); ++ivar )
  {
    MnvH1D *xsec  = utils->GetXSecHist( fRawXSec, targetZ, *ivar, helicity, targetID );
    xsecHists.push_back( xsec );
  }

     for( vector<string>::iterator ivar = vars.begin(); ivar != vars.end(); ++ivar ){
    

       const string& var = *ivar;
       f1->cd();     


       MnvH1D *dsignal = (MnvH1D*)f1->Get(Form("data_%s", var.c_str()));
        if( !dsignal )
            cout << "ERROR: could not get Data MnvH1D for variable: " << *ivar << endl;
        dHists.push_back( dsignal );
   
        MnvH1D *dunfold = (MnvH1D*)f1->Get(Form("unfolded_data_%s", var.c_str()));
        if( !dunfold )
            cout << "ERROR: could not get Unfolded Data MnvH1D for variable: " << *ivar << endl;
        dHists_unfolded.push_back( dunfold );
        
        MnvH1D *mcsignal = (MnvH1D*)f1->Get(Form("signal_%s", var.c_str()));
        if( !mcsignal )
            cout << "ERROR: could not get MC MnvH1D for variable: " << *ivar << endl;
        mcHists.push_back( mcsignal );

        MnvH1D *mcunfold = (MnvH1D*)f1->Get(Form("unfolded_signal_%s", var.c_str()));
        if( !mcunfold )
            cout << "ERROR: could not get Unfolded MC MnvH1D for variable: " << *ivar << endl;
        mcHists_unfolded.push_back( mcunfold );
        
        MnvH1D *effNum = (MnvH1D*)Eff->Get(Form("h_mc_%s", var.c_str()));
        MnvH1D *effDen = (MnvH1D*)Eff->Get(Form("h_truth_%s", var.c_str()));
        effNum->Divide(effNum, effDen);
        effHists.push_back( effNum );
    }  


    double sumScaleFactor = 1.0;
    for( unsigned int i = 0; i != vars.size(); ++i )
    {
        const string var = vars[i];

    if( !dHists_unfolded[i] || !mcHists_unfolded[i] || !effHists[i] )
    {
      cout << "  Not enough histograms for variable " << var << endl;
      continue;
    }

    MnvH1D *mnvDunfolded = dynamic_cast<MnvH1D*>( dHists_unfolded[i]->Clone( "tmpDUnfold" ) );
    MnvH1D *mnvMCunfolded = dynamic_cast<MnvH1D*>( mcHists_unfolded[i]->Clone( "tmpMCUnfold" ) );
    mnvDunfolded->ClearSysErrorMatrices();
    mnvMCunfolded->ClearSysErrorMatrices();
      
      cout<<i<<"  185"<<endl;

    MnvH1D *mnvEff = dynamic_cast<MnvH1D*>( effHists[i]->Clone( "tmpEff" ) );


    //scale all by N nucleons
    mnvDunfolded->Scale( 1. / nnuc_D );
    mnvMCunfolded->Scale( 1. / nnuc_MC );

    //do efficiency correction on unfolded distributions
    MnvH1D *mnvDEffCorr = utils->EfficiencyCorrect( mnvDunfolded, mnvEff );
    MnvH1D * mnvMCEffCorr = utils->EfficiencyCorrect( mnvMCunfolded, mnvEff );

   
    //divide by flux to get cross section
    string tmpplaylist= playlist;
    if(tmpplaylist=="AllNuME") tmpplaylist="minervame1A";
    MnvH1D *mnvDXSec = (MnvH1D*)mnvDEffCorr->Clone("tmpXSec_D");
    if( "Enu" == var )
      utils->DivideByFlux( mnvDXSec, tmpplaylist );
    else
      utils->DivideByIntegratedFlux( mnvDXSec, tmpplaylist );

    MnvH1D *mnvMCXSec = (MnvH1D*)mnvMCEffCorr->Clone("tmpXSec_MC");
    if( "Enu" == var )
      utils->DivideByFlux( mnvMCXSec, tmpplaylist );
    else
      utils->DivideByIntegratedFlux( mnvMCXSec, tmpplaylist );
 
      MnvH1D *generatedXSec = xsecHists[i];
//      MnvH1D *generatedEvtHist = evtrateHists[i];

    mnvDXSec->Scale( 1.0E38 / datapot );
    mnvMCXSec->Scale( 1.0E38 / datapot);
    mnvDEffCorr->Scale( 1.0E38 / datapot );
    mnvMCEffCorr->Scale( 1.0E38 / datapot );

    //generatedXSec->Scale(1.0E38);
    cout<<i<<"  220"<<endl;

    
    const double minVal = binsDef->GetDISVarMinVal( var );
    const double maxVal = binsDef->GetDISVarMaxVal( var );

    //const double minVal = NukeUtils::Get().GetVarMinVal( var );
    //const double maxVal = NukeUtils::Get().GetVarMaxVal( var );
    mnvDunfolded->GetXaxis()->SetRangeUser( minVal, maxVal );
    mnvMCunfolded->GetXaxis()->SetRangeUser( minVal, maxVal );
    mnvDEffCorr->GetXaxis()->SetRangeUser( minVal, maxVal );
    mnvMCEffCorr->GetXaxis()->SetRangeUser( minVal, maxVal );
    mnvDXSec->GetXaxis()->SetRangeUser( minVal, maxVal );
    mnvMCXSec->GetXaxis()->SetRangeUser( minVal, maxVal );
    if( generatedXSec )
    {
      generatedXSec->GetXaxis()->SetRangeUser( minVal, maxVal );
      generatedXSec->SetMarkerColor( kBlue-2 );
      generatedXSec->SetLineColor( kBlue-2);
      generatedXSec->SetMarkerStyle( 34 );
      generatedXSec->SetMarkerSize( 2.5 );
    }

    cout << "   Making dataMC plots of " << var << " for (targetID, targetZ, helicity) = "
      << targetID << ", " << targetZ << " = " << targetString << ", " << playlist << endl;

    const string x_title = binsDef->GetXaxisTitle( var );
    const string y_title = binsDef->GetYaxisTitle( var );
    const std::string xsec_label = utils->GetTotalXSecString(targetID, targetZ) + utils->GetTotalXSecUnits();

    double potX = .28, potY = .83;
    double potXLog = .68, potYLog = .76;
    double titleSize = .038;
    if( targetID > 10 )
      titleSize = .032;
   

        cout << "after unfolding now how many unfolded? type? " << mnvDunfolded->Integral() << ", " << FileType::kAny << endl;
        cout << "after unfolding now how many unfoldedData? type? " << mnvMCunfolded->Integral() << ", " << FileType::kAny << endl;
        for( int bin=1; bin < mnvDunfolded->GetNbinsX()+1; bin++ ){
        cout << "bin content, bin error? type? " << bin << ", " << mnvDunfolded->GetBinContent(bin) << ", " << mnvDunfolded->GetBinError(bin) << ", " << FileType::kAny << endl;
        cout << "bin content, bin error? type? " << bin << ", " << mnvMCunfolded->GetBinContent(bin) << ", " << mnvMCunfolded->GetBinError(bin) << ", " << FileType::kAny << endl;
}
    //====== Draw-MC Unfolded Signal ======//
    {
 
      TString cName = Form("EnergyPlotUnfoldedPerNucleon_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas c( cName, cName, 1280, 800 );
      mnvMCunfolded->GetXaxis()->SetTitle( x_title.c_str() );
      mnvMCunfolded->GetYaxis()->SetTitle( string( y_title + " / nucleon").c_str() );
      cout<<"can you reach here"<<endl;     
      //mnvPlotter.DrawDataMCWithErrorBand( mnvDunfolded, mnvMCunfolded, dataMCScale, "TC-TR");
      mnvPlotter.DrawDataMCWithErrorBand( mnvDunfolded, mnvMCunfolded, dataMCScale, "TC-TR", false, NULL, NULL, false, true);
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary("TL");
      if( WRITE_POT ) mnvPlotter.WriteNorm("POT-Normalized", potX, potY, .04, datapot );
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Unfolded - %s", targetString.c_str() ), titleSize );
      mnvPlotter.MultiPrint( &c );

      cName = Form( "ResultPlotUnfoldedRatioPerNucleon_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas cR( cName, cName, 1280, 800 );
      //mnvMCunfolded->GetYaxis()->SetTitle( "Data/MC" );
      TString Yaxis = "Data/MC"; 
      std::string x_label = var.c_str();
      mnvPlotter.DrawDataMCRatio( mnvDunfolded, mnvMCunfolded, dataMCScale, false, true, 0.5, 1.5, x_label.c_str(), Yaxis );
      mnvPlotter.AddChi2Label( mnvDunfolded, mnvMCunfolded, dataMCScale, "TC", .06);
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary("TR");
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Unfolded - %s", targetString.c_str() ), titleSize );
      mnvPlotter.MultiPrint( &cR );
    }
    //if( NukeUtils::Get().ShowVarAsLog( var ) )
    {
      TString cName = Form("ResultPlotUnfoldedLogPerNucleon_%s_%s", var.c_str(), suffix.c_str()  );
      TCanvas c( cName, cName, 1280, 800 );
      c.cd(0)->SetLogy(1);
      mnvMCunfolded->GetXaxis()->SetTitle( x_title.c_str() );
      mnvMCunfolded->GetYaxis()->SetTitle( string( y_title + " / nucleon").c_str() );
      //mnvPlotter.DrawDataMCWithErrorBand( mnvDunfolded, mnvMCunfolded, dataMCScale, "TC-TR");
      mnvPlotter.DrawDataMCWithErrorBand( mnvDunfolded, mnvMCunfolded, dataMCScale, "TC-TR", false, NULL, NULL, false, true);
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary("BL");
      if( WRITE_POT ) mnvPlotter.WriteNorm("POT-Normalized",  potXLog, potYLog, .04, datapot );
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Unfolded - %s", targetString.c_str() ), titleSize );
      mnvPlotter.MultiPrint( &c );
    }
    //========== End of Data-MC Unfolded ============//



    //====== Cross Section ======//
    {
      TString cName = Form( "ResultPlotCrossSection_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas c( cName, cName, 1280, 800 );
      mnvMCXSec->GetXaxis()->SetTitle( x_title.c_str() );
      mnvMCXSec->GetYaxis()->SetTitle( xsec_label.c_str() );
      //mnvPlotter.DrawDataMCWithErrorBand( mnvDXSec, mnvMCXSec, dataMCScale, "TC-TR");
      mnvPlotter.DrawDataMCWithErrorBand( mnvDXSec, mnvMCXSec, dataMCScale, "TC-TR", false, NULL, NULL, false, true);
      mnvPlotter.AddChi2Label( mnvDXSec, mnvMCXSec, dataMCScale, "TL", .06);
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary("TL");
      if( WRITE_POT ) mnvPlotter.WriteNorm("",  potX, potY, .04, datapot );
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Cross Section - %s", targetString.c_str() ), titleSize );
      mnvPlotter.MultiPrint( &c );

      cName = Form("ResultPlotCrossSectionRatio_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas cR( cName, cName, 1280, 800 );
      //mnvMCXSec->GetYaxis()->SetTitle( "Data/MC" );
      TString Yaxis = "Data/MC"; 
      std::string x_label = var.c_str();
      mnvPlotter.DrawDataMCRatio(  mnvDXSec, mnvMCXSec, dataMCScale, false, true, 0.0, 2.0, x_label.c_str(), Yaxis );
      mnvPlotter.AddChi2Label( mnvDXSec, mnvMCXSec, dataMCScale, "TC",.06);
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary("TR");
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Cross Section - %s", targetString.c_str() ), titleSize );
      mnvPlotter.MultiPrint( &cR );
    }
    //if( NukeUtils::Get().ShowVarAsLog( var ) )
    {
      TString cName = Form("ResultPlotCrossSectionLog_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas c( cName, cName, 1280, 800 );
      c.cd(0)->SetLogy(1);
      mnvMCXSec->GetXaxis()->SetTitle( x_title.c_str() );
      mnvMCXSec->GetYaxis()->SetTitle(xsec_label.c_str() );
      //mnvPlotter.DrawDataMCWithErrorBand( mnvDXSec, mnvMCXSec, dataMCScale, "TC-TR");
      mnvPlotter.DrawDataMCWithErrorBand( mnvDXSec, mnvMCXSec, dataMCScale, "TC-TR", false, NULL, NULL, false, true);
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary("BL");
      if( WRITE_POT ) mnvPlotter.WriteNorm("",  potXLog, potYLog, .04, datapot );
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Cross Section - %s", targetString.c_str() ), titleSize );
      mnvPlotter.MultiPrint( &c );
    }


    //Cross Section to True
    if( generatedXSec )
    {
      const int old_dataMarker = mnvPlotter.data_marker;
      const int old_dataColor  = mnvPlotter.data_color;
      const double old_dataSize = mnvPlotter.data_marker_size;

      mnvPlotter.data_marker = 34;
      mnvPlotter.data_color  = kBlue-2;
      mnvPlotter.data_marker_size = 2.5;

      TString cName = Form("ResultPlotMCGenerated_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas c( cName, cName, 1280, 800 );
      mnvMCXSec->SetTitle("MC");
      generatedXSec->SetTitle("GENIE Calculation");
      mnvMCXSec->GetXaxis()->SetTitle( x_title.c_str() );
      mnvMCXSec->GetYaxis()->SetTitle( xsec_label.c_str() );
      //mnvPlotter.DrawDataMCWithErrorBand( generatedXSec, mnvMCXSec, dataMCScale, "TC-TR", true);
      mnvPlotter.DrawDataMCWithErrorBand( generatedXSec, mnvMCXSec, dataMCScale, "TC-TR", true, NULL, NULL, false, true);
      mnvPlotter.AddChi2Label( generatedXSec, mnvMCXSec, dataMCScale, "TL", .06);
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary("TL");

//      generatedXSec->SetLineWidth(4);
//      generatedXSec->SetLineColor(kBlue-2);
//      generatedXSec->Draw("SAME");

      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Cross Section - %s", targetString.c_str() ), titleSize );
      mnvPlotter.MultiPrint( &c );

      cName = Form("ResultPlotMCGeneratedRatio_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas cR( cName, cName, 1280, 800 );
      generatedXSec->GetYaxis()->SetTitle( "MC / GENIE Calculation" );
      //mnvPlotter.DrawDataMCRatio(  mnvMCXSec, generatedXSec, 1./dataMCScale, true, false, true, -1., -1., "MC / GENIE Calculation" );//not sure thats really what I want 
      mnvPlotter.DrawDataMCRatio(  mnvMCXSec, generatedXSec, 1./dataMCScale, true, false, -1., -1., "MC / GENIE Calculation" );//not sure thats really what I want 
      mnvPlotter.AddChi2Label( mnvMCXSec, generatedXSec,  1./dataMCScale, "TC", .06);
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary("TR");
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Cross Section - %s", targetString.c_str() ), titleSize );
      mnvPlotter.MultiPrint( &cR );
        

      //return data marker style
      mnvPlotter.data_marker = old_dataMarker;
      mnvPlotter.data_color  = old_dataColor;
      mnvPlotter.data_marker_size = old_dataSize;
    }


    //===============================
    // ERRORS
    //===============================
    mnvPlotter.legend_n_columns = 2;
    mnvPlotter.height_nspaces_per_hist = 1.;

    string prelimLoc = "TR";
    string legLoc    = "TL";

    //=======================
    // Unfolded
    //Draw an error summary for unfolded signal
    {
      TString cName = Form("ResultPlotUnfoldedErrSummary_Data_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas cErrSummary( cName, cName, 1280,920 );
      //mnvMCunfolded->GetXaxis()->SetTitle( x_title.c_str() );
      mnvMCunfolded->GetYaxis()->SetTitle( "Fractional Uncertainty" );
      std::string x_label = var.c_str(); 
      PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
      plotter.DrawErrorSummary( mnvDunfolded, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
      //plotter.DrawErrorSummary( mnvDunfolded, legLoc, false );
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary( prelimLoc, .04 );
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Errors on Unfolded Signal in %s", targetString.c_str() ), titleSize);
      mnvPlotter.MultiPrint( &cErrSummary );
    }

    //Draw an error summary for unfolded signal.  no stat errors
    {
      TString cName = Form("ResultPlotUnfoldedSysErrSummary_MC_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas cErrSummary( cName, cName, 1280,920 );
      mnvMCunfolded->GetXaxis()->SetTitle( x_title.c_str() );
      mnvMCunfolded->GetYaxis()->SetTitle( "Fractional Uncertainty" );
     std::string x_label = var.c_str(); 
      PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
      plotter.DrawErrorSummary( mnvMCunfolded, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
      //plotter.DrawErrorSummary( mnvMCunfolded, legLoc, false );
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary( prelimLoc, .04 );
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Sys. Errors on Unfolded Signal in %s", targetString.c_str() ), titleSize);
      mnvPlotter.MultiPrint( &cErrSummary );
    }


    //=======================
    // Eff. Corrected
    // Draw an error summary for eff. corrected signal
    {
      TString cName = Form("ResultPlotEffCorrErrSummary_Data_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas cErrSummary( cName, cName, 1280,920 );
      mnvMCEffCorr->GetXaxis()->SetTitle( x_title.c_str() );
      mnvMCEffCorr->GetYaxis()->SetTitle( "Fractional Uncertainty" );
     std::string x_label = var.c_str(); 
      PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
      plotter.DrawErrorSummary( mnvDEffCorr, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
      //plotter.DrawErrorSummary(mnvDEffCorr, legLoc, false );
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary( prelimLoc, .04 );
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Errors on Eff. Corrected Data Signal in %s", targetString.c_str() ), titleSize);
      mnvPlotter.MultiPrint( &cErrSummary );
    }

    //Draw an error summary for unfolded signal.  no stat errors
    {
      TString cName = Form("ResultPlotEffCorrSysErrSummary_MC_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas cErrSummary( cName, cName, 1280,920 );
      mnvMCEffCorr->GetXaxis()->SetTitle( x_title.c_str() );
      mnvMCEffCorr->GetYaxis()->SetTitle( "Fractional Uncertainty" );
      PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
      plotter.DrawErrorSummary( mnvMCEffCorr, legLoc, false );
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary( prelimLoc, .04 );
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Sys. Errors on Eff. Corrected Unfolded MCSignal in %s", targetString.c_str() ), titleSize);
      mnvPlotter.MultiPrint( &cErrSummary );
    }


    //=======================
    // Cross Section
    {
      TString cName = Form("ResultPlotCrossSectionErrSummary_Data_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas cErrSummary( cName, cName, 1280,920 );
      mnvMCEffCorr->GetXaxis()->SetTitle( x_title.c_str() );
      mnvMCEffCorr->GetYaxis()->SetTitle( "Fractional Uncertainty" );
      PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
      plotter.DrawErrorSummary(  mnvDEffCorr, legLoc, false );
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary( prelimLoc, .04 );
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Errors on Data Cross Section in  %s", targetString.c_str() ), titleSize);
      mnvPlotter.MultiPrint( &cErrSummary );
    }

    {
      TString cName = Form("ResultPlotCrossSectionSysErrSummary_MC_%s_%s", var.c_str(), suffix.c_str() );
      TCanvas cErrSummary( cName, cName, 1280,920 );
      mnvMCEffCorr->GetXaxis()->SetTitle( x_title.c_str() );
      mnvMCEffCorr->GetYaxis()->SetTitle( "Fractional Uncertainty" );
      PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
      plotter.DrawErrorSummary( mnvMCEffCorr, legLoc, false );
      if( WRITE_PRELIMINARY )
        mnvPlotter.WritePreliminary( prelimLoc, .04 );
      if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Sys. Errors on MC Cross Section in %s", targetString.c_str() ), titleSize);
      mnvPlotter.MultiPrint( &cErrSummary );
    }

    mnvPlotter.height_nspaces_per_hist = 1.2;
    mnvPlotter.legend_n_columns = 1;


    //cleanup
    cout << "   ...cleaning up..." << endl;
    delete mnvDEffCorr;
    delete mnvMCEffCorr;
    delete mnvDunfolded;
    delete mnvMCunfolded;
    delete mnvEff;
    delete generatedXSec;
    mnvPlotter.CleanTmp();
  }

  //cleanup
  delete f1;
  delete fRawXSec;


  return 0;



//}//variable loop




}// End of Main





void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name ){
  
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


