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
       effFileName = Form("%s/Hists_Efficiency_t%d_z%02d_Nu_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str() ); 
       eventfile = Form("%s/Hists_EventSelection_t%d_z%02d_Nu_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str() ); 
     }
   else{
       histFileName = Form("%s/Unfolded_t%d_z%02d_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str() ); 
       effFileName = Form("%s/Hists_Efficiency_t%d_z%02d_Nu_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str() ); 
       eventfile = Form("%s/Hists_EventSelection_t%d_z%02d_Nu_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str() ); 
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
   vars.push_back("Emu");
   vars.push_back("Enu");
//   vars.push_back("Ehad");
//   vars.push_back("x");
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
 
   
   TParameter<double> *mcPOT = (TParameter<double>*)f0->Get("MCPOT");
   TParameter<double> *dataPOT = (TParameter<double>*)f0->Get("DataPOT");
   double mcpot = mcPOT->GetVal();
   double datapot = dataPOT->GetVal(); 
   double dataMCScale = datapot/mcpot;
   cout<<"MCPOT = "<<mcpot<<"DataPOT = "<< datapot << "Scale = " << dataMCScale<<endl;
   
   TObjArray hlistUntunedTrans, hlistUntunedContin, hlisttunedTrans, hlisttunedContin;


   std::vector<MnvH1D*> dBGHists, mcBGHists, dHistsBS, mcHistsBS, mcHistsWrongTarget, mcHistsNonDIS;
   std::vector<MnvH1D*> dHists, dHists_unfolded;
   std::vector<MnvH1D*> mcHists, mcHists_unfolded;
   vector< MnvH1D*> hists_Trans_data, hists_Trans_in_Trans, hists_Contin_in_Trans, hists_Signal_in_Trans;
   vector< MnvH1D*> hists_Contin_data, hists_Trans_in_Contin, hists_Contin_in_Contin, hists_Signal_in_Contin;
   vector< MnvH1D*> histos_mc_trans, histos_mc_contin, data_sideband, data_sideband_Contin;
   std::vector<MnvH1D*> effHists;


    const bool useBG = targetID < 10;

     for( vector<string>::iterator ivar = vars.begin(); ivar != vars.end(); ++ivar ){
    
       const string& var = *ivar;
       f1->cd();     
       f0->cd();
     //Save the before background subtracted distributions

       MnvH1D *dataBS = (MnvH1D*)f0->Get(Form("h_data_%s", var.c_str()));
        if( !dataBS )
            cout << "ERROR: could not get Data Before Subtraction MnvH1D for variable: " << *ivar << endl;
        dHistsBS.push_back( dataBS );
  
        MnvH1D *mcBS = (MnvH1D*)f0->Get(Form("h_mc_%s", var.c_str()));
        if( !mcBS )
            cout << "ERROR: could not get MC before subtraction MnvH1D for variable: " << *ivar << endl;
        mcHistsBS.push_back( mcBS );

     
       MnvH1D *BGData = (MnvH1D*)f1->Get(Form("BG_data_%s", var.c_str()));
        if( !BGData )
            cout << "ERROR: could not get Data Before Subtraction MnvH1D for variable: " << *ivar << endl;
        dBGHists.push_back( BGData );
  
        MnvH1D *BGMC = (MnvH1D*)f1->Get(Form("BG_mc_%s", var.c_str()));
        if( !BGMC )
            cout << "ERROR: could not get MC before subtraction MnvH1D for variable: " << *ivar << endl;
        mcBGHists.push_back( BGMC );

        MnvH1D *MCWrongTarget = (MnvH1D*)f1->Get(Form("BG_WrongTarget_%s", var.c_str()));
        if( !MCWrongTarget )
            cout << "ERROR: could not get MC before subtraction MnvH1D for variable: " << *ivar << endl;
        mcHistsWrongTarget.push_back( MCWrongTarget );

        MnvH1D *MCNonDIS = (MnvH1D*)f1->Get(Form("BG_NonDIS_%s", var.c_str()));
        if( !MCNonDIS )
            cout << "ERROR: could not get MC before subtraction MnvH1D for variable: " << *ivar << endl;
        mcHistsNonDIS.push_back( MCNonDIS );
      ///////////////////////////////////////////////////


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
        
        cout<<"MC Integral   "<<mcHists[i]->Integral()<<"  "<<mcHists_unfolded[i]->Integral()<<endl;
        cout<<"Data Integral "<<dHists[i]->Integral()<<"  "<<dHists_unfolded[i]->Integral()<<endl;
        
        cout << "   Making dataMC plots of " << var << " for (targetID, targetZ, helicity) = "<< targetID << ", " << targetZ << " = " << endl;
        
        if( ! ( dHists[i] && mcHists[i] ) )
        {
            cout << "  Not enough histograms for variable " << var << endl;
            continue;
        }
        
        
        //Get the Sum of the energy hists
        int sumZ = -1;
        
        if(useBG)
            sumZ = targetZ;
        MnvH1D *sumMCSig; 
        MnvH1D *sumDataSig;
        
        MnvH1D *sumMCSigUF;
        MnvH1D *sumDataSigUF;   
        
        //We've already corrected for the non-DIS contam at this point
       /* if(sumMaterial){
            sumMCSig = SumMaterialHists("Unfolded", "signal", var, mc_type, sumZ, false);
            sumDataSig = SumMaterialHists("Unfolded", "signal", var, data_type, sumZ, false);
            
            sumMCSigUF = SumMaterialHists("Unfolded", "unfolded", var, mc_type, sumZ, false);
            sumDataSigUF = SumMaterialHists("Unfolded", "unfolded", var, data_type, sumZ, false);
        }
        */
/*        if(var == "x" && sumMaterial){
            const vector<double> xBins = binsDef->GetEnergyBins("x");
            const unsigned int nBins = xBins.size();
            
            double Events = 0;
            for(unsigned int bin = 0; bin != xBins.size()-1; ++bin){
                //Events[bin] += dataMCScale*carbonHist[0]->GetBinContent(bin+1);
                
                mcxBjTable <<  setprecision(2) << xBins[bin] << " - " << xBins[bin+1] << setprecision(0) << std::fixed << " & " << sumDataSigUF->GetBinContent(bin+1) << endl;
                Events += sumDataSigUF->GetBinContent(bin+1);
                
            }
            mcxBjTable <<  "TOTAL: & " << setprecision(0) << Events << endl; 
            
            
        }
        
  */      
        
        const double minValDIS = binsDef->GetDISVarMinVal( var );
        const double maxValDIS = binsDef->GetDISVarMaxVal( var );
        
        ///const double minVal = NukeUtils::Get().GetVarMinVal( var );
        ///const double maxVal = NukeUtils::Get().GetVarMaxVal( var );
        
        MnvH1D *mnvBGD  = dBGHists[i];
        MnvH1D *mnvBGMC  = mcBGHists[i];
    
        //Background for breakdown
        MnvH1D *mnvWrongTarget = mcHistsWrongTarget[i];
        MnvH1D *mnvNonDIS = mcHistsNonDIS[i];
         
        MnvH1D *mnvDBS  = dHistsBS[i];
        MnvH1D *mnvMCBS  = mcHistsBS[i];
        
        MnvH1D *mnvD  = dHists[i];
        MnvH1D *mnvMC = mcHists[i];
        
        MnvH1D *mnvD_uf  = dHists_unfolded[i];
        MnvH1D *mnvMC_uf = mcHists_unfolded[i];
        
        MnvH1D *mnvD_unf_effcorr  = dHists_unfolded[i];
        MnvH1D *mnvMC_unf_effcorr = mcHists_unfolded[i];
        
        MnvH1D *eff = effHists[i];
        
        for( int bin=0; bin < mnvBGD->GetNbinsX(); bin++ ){
        cout << "bin content, bin error? type? " << bin << ", " << mnvBGD->GetBinContent(bin) << ", " << mnvBGD->GetBinError(bin) << ", " << FileType::kAny << endl;
        cout << "bin content, bin error? type? " << bin << ", " << mnvD->GetBinContent(bin) << ", " << mnvD->GetBinError(bin) << ", " << FileType::kAny << endl;
        }
        
        cout<<"Spare my life plz for MC"<<endl;
        for( int bin=0; bin < mnvBGMC->GetNbinsX(); bin++ ){
        cout << "bin content, bin error? type? " << bin << ", " << mnvBGMC->GetBinContent(bin) << ", " << mnvBGMC->GetBinError(bin) << ", " << FileType::kAny << endl;
        cout << "bin content, bin error? type? " << bin << ", " << mnvMC->GetBinContent(bin) << ", " << mnvMC->GetBinError(bin) << ", " << FileType::kAny << endl;
        }
        
        mnvD_unf_effcorr = utils->EfficiencyCorrect(mnvD_unf_effcorr,eff);
        mnvMC_unf_effcorr = utils->EfficiencyCorrect(mnvMC_unf_effcorr,eff);
        cout<<"Bin    Uncorr  Eff   Ratio   Corr "<<endl;
        for (int i=0; i<mnvD_unf_effcorr->GetXaxis()->GetNbins(); i++){
            
            cout<<i<<"  "<<mnvD_uf->GetBinContent(i)<<"  "<<eff->GetBinContent(i)<<"  "<<mnvD_uf->GetBinContent(i)/eff->GetBinContent(i)<<"   "<<mnvD_unf_effcorr->GetBinContent(i)<<endl;
            
        }
        
      
            double Nucleons = -999.;
            
            if(targetZ==6)
                Nucleons =  TargetUtils::Get().GetPassiveTargetNNucleons( 3, targetZ, 0 );
            if(targetZ==26)
                Nucleons =  TargetUtils::Get().GetPassiveTargetNNucleons( targetID, targetZ, 0 );
            if(targetZ==82)
                Nucleons =  TargetUtils::Get().GetPassiveTargetNNucleons( targetID, targetZ, 0 );
            if(targetID>10 )
                Nucleons =  TargetUtils::Get().GetTrackerNNucleons( 12.0, 0);
       
        
        MnvH1D* mnvD_unf_effcorr_nucleons = (MnvH1D*) mnvD_unf_effcorr->Clone("mnvD_unf_effcorr_nucleons");
        MnvH1D* mnvMC_unf_effcorr_nucleons = (MnvH1D*) mnvMC_unf_effcorr->Clone("mnvMC_unf_effcorr_nucleons");
        
        
        mnvD_unf_effcorr_nucleons->Scale(1./Nucleons);
        mnvMC_unf_effcorr_nucleons->Scale(1./Nucleons);

        
        mnvD->GetXaxis()->SetRangeUser( minValDIS, maxValDIS );
        mnvMC->GetXaxis()->SetRangeUser( minValDIS, maxValDIS );
        mnvD_uf->GetXaxis()->SetRangeUser( minValDIS, maxValDIS );
        mnvMC_uf->GetXaxis()->SetRangeUser( minValDIS, maxValDIS );
       
        const string x_titleTrue = "Unfolded " + binsDef->GetXaxisTitle( var );
//        const string x_titleTrue = NukeUtils::Get().GetXaxisTitle( var );
        const string x_titleReco = "Reconstructed " + binsDef->GetXaxisTitle( var );
        const string y_title = binsDef->GetYaxisTitle( var );
        
        string legPos = "TR";
        string prelimPos = "TL", prelimPosR = "TL", prelimPosLog = "BL";
        double potX = .3, potY = .85;
        double titleSize = .04;
        if( targetID > 10 )
            titleSize = .035;
        
        if( var == "PhiMu" )
        {
            potX = .3;
            potY = .8;
        }
   
        mnvPlotter.height_nspaces_per_hist = 1.5;
        mnvPlotter.legend_n_columns = 1;
        //gStyle->SetEndErrorSize(0.0);
        if( GRAYSCALE ){
            mnvPlotter.mc_error_style = 1001;
            mnvPlotter.data_line_width = 4;
        }
        
    /*    //find the maximum y-value for comparing plots
        if(sumMaterial){
            sumMCSig->GetXaxis()->SetRangeUser(minValDIS, maxValDIS);
            sumMCSigUF->GetXaxis()->SetRangeUser(minValDIS, maxValDIS);
        }
        double maxY = 0.0;
        //The fudge factors roughly do the bin width normalizaion, and give some headroom for labels, etc. 
   cout<<"Can i see here"<<endl;            
        if(sumMaterial){
            if( var != "x" )
                maxY = 0.33*max( dataMCScale*sumMCSig->GetBinContent(sumMCSig->GetMaximumBin()), dataMCScale*sumMCSigUF->GetBinContent(dataMCScale*sumMCSigUF->GetMaximumBin()) );
            else
                maxY = 2000;
//                maxY = 10.0*max( dataMCScale*sumMCSig->GetBinContent(sumMCSig->GetMaximumBin()), dataMCScale*sumMCSigUF->GetBinContent(dataMCScale*sumMCSigUF->GetMaximumBin()) );
            
        }*/
        //maxY = 2000;
        cout<<"269"<<endl;
        cout<<"MC Integral   "<<mnvMC->Integral()<<"  "<<mnvMC_uf->Integral()<<endl;
        cout<<"Data Integral "<<mnvD->Integral()<<"  "<<mnvD_uf->Integral()<<endl;
        
        

//}
//////////////For before background plots /////////////////////////////
/*
        {
	  TCanvas ctest;
	  MnvH1D* sig2=mnvMCBS->Clone();
	  cout<<"sig2 nombinwidth "<<sig2->GetNormBinWidth();
	  sig2->Scale(sig2->GetNormBinWidth(), "width");
	  sig2->SaveAs(Form("signal_%s.png", var.c_str() ));
            TString cName = Form("DIS_BeforeBG_breakdown_dataMC_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cDIS( cName, cName, 1280, 800 );
            cDIS.SetGrayscale( GRAYSCALE );
            
            mnvMCBS->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvMCBS->GetYaxis()->SetTitle( y_title.c_str() );
            mnvMCBS->GetXaxis()->SetRangeUser(minValDIS, maxValDIS);
            cout<<"280"<<endl;
            mnvMCBS->SetMaximum(2000.0);
            mnvDBS->SetMaximum(2000.0);
            

            if(var == "x"){
                mnvDBS->SetNormBinWidth(0.1);
                mnvMCBS->SetNormBinWidth(0.1);
                mnvBGD->SetNormBinWidth(0.1);
                mnvBGMC->SetNormBinWidth(0.1);
            }

            if(var == "Enu"){
                mnvDBS->SetNormBinWidth(1.0);
                mnvMCBS->SetNormBinWidth(1.0);
                mnvBGD->SetNormBinWidth(1.0);
                mnvBGMC->SetNormBinWidth(1.0);
            }

            if(var == "y"){
                mnvDBS->SetNormBinWidth(0.1);
                mnvMCBS->SetNormBinWidth(0.1);
                mnvBGD->SetNormBinWidth(0.1);
                mnvBGMC->SetNormBinWidth(0.1);
            }

	    mnvPlotter.DrawDataMCWithErrorBand( mnvDBS, mnvMCBS, dataMCScale, legPos, false, mnvWrongTarget, mnvNonDIS, mnvBGD, false, true);
 
            if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPos);
            if( WRITE_POT ) mnvPlotter.WriteNorm("POT-Normalized", potX, potY, .04, datapot );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "DIS Signal - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cDIS );
        }*/
        {
	  TCanvas ctest;
	  MnvH1D* sig2=mnvMCBS->Clone();
	  cout<<"sig2 nombinwidth "<<sig2->GetNormBinWidth();
	  sig2->Scale(sig2->GetNormBinWidth(), "width");
	  sig2->SaveAs(Form("signal_%s.png", var.c_str() ));
            TString cName = Form("DIS_BeforeBG_dataMC_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cDIS( cName, cName, 1280, 800 );
            cDIS.SetGrayscale( GRAYSCALE );
            
            mnvMCBS->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvMCBS->GetYaxis()->SetTitle( y_title.c_str() );
            mnvMCBS->GetXaxis()->SetRangeUser(minValDIS, maxValDIS);
            cout<<"280"<<endl;
            mnvMCBS->SetMaximum(2000.0);
            mnvDBS->SetMaximum(2000.0);
            

            if(var == "x"){
                mnvDBS->SetNormBinWidth(0.1);
                mnvMCBS->SetNormBinWidth(0.1);
                mnvBGD->SetNormBinWidth(0.1);
                mnvBGMC->SetNormBinWidth(0.1);
            }

            if(var == "Enu"){
                mnvDBS->SetNormBinWidth(1.0);
                mnvMCBS->SetNormBinWidth(1.0);
                mnvBGD->SetNormBinWidth(1.0);
                mnvBGMC->SetNormBinWidth(1.0);
            }

            if(var == "y"){
                mnvDBS->SetNormBinWidth(0.1);
                mnvMCBS->SetNormBinWidth(0.1);
                mnvBGD->SetNormBinWidth(0.1);
                mnvBGMC->SetNormBinWidth(0.1);
            }


	    mnvPlotter.DrawDataMCWithErrorBand( mnvDBS, mnvMCBS, dataMCScale, legPos, false, mnvBGMC, mnvBGD, false, true);
	    //mnvPlotter.DrawDataMCWithErrorBand( mnvDBS, mnvMCBS, dataMCScale, legPos, false, NULL, NULL, false, true);
	    //mnvPlotter.DrawDataMCWithErrorBand( mnvDBS, mnvMCBS, dataMCScale, legPos, false, mnvNonDIS, mnvBGD, false, true);
 
            if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPos);
            if( WRITE_POT ) mnvPlotter.WriteNorm("POT-Normalized", potX, potY, .04, datapot );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "DIS Signal - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cDIS );
        }

        {
            
            TString cName = Form( "DIS_BeforeBG_dataMCRatio_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cR( cName, cName, 1280, 800 );      
            cR.SetGrayscale( GRAYSCALE );
            mnvMCBS->GetYaxis()->SetTitle( "Data/MC" );
            TString Yaxis = "Data/MC"; 
             std::string x_label = var.c_str();
	    //mnvPlotter.DrawDataMCRatio( mnvD, mnvMC, dataMCScale, false, true, true, .6, 1.7 );
	    mnvPlotter.DrawDataMCRatio( mnvDBS, mnvMCBS, dataMCScale, true, true, 0.0, 2.0, x_label.c_str(), Yaxis );
            mnvPlotter.AddChi2Label(mnvDBS, mnvMCBS, dataMCScale, "TL", .06, 0.0, true );
	    if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPosR);
            if( WRITE_POT ) 
                mnvPlotter.WriteNorm("POT-Normalized", "TR", .04, datapot );
            else
                mnvPlotter.WriteNorm("Shape Only", "TR", .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Before BG Subtraction Signal - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cR );
            
        }





//////////////////////////////////////////////////////////////////////





        //====== Draw Data-MC Overlay Signal ======//
        {
	  TCanvas ctest;
	  MnvH1D* sig2=mnvMC->Clone();
	  cout<<"sig2 nombinwidth "<<sig2->GetNormBinWidth();
	  sig2->Scale(sig2->GetNormBinWidth(), "width");
	  sig2->SaveAs(Form("signal_%s.png", var.c_str() ));
            TString cName = Form("DIS_SignalPlot_dataMC_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cDIS( cName, cName, 1280, 800 );
            cDIS.SetGrayscale( GRAYSCALE );
            
            mnvMC->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvMC->GetYaxis()->SetTitle( y_title.c_str() );
            mnvMC->GetXaxis()->SetRangeUser(minValDIS, maxValDIS);
            cout<<"280"<<endl;
            mnvMC->SetMaximum(2000.0);
            mnvD->SetMaximum(2000.0);
            
	    mnvPlotter.DrawDataMCWithErrorBand( mnvD, mnvMC, dataMCScale, legPos, false, NULL, NULL, false, true);
	    //mnvPlotter.DrawDataMCWithErrorBand( mnvD, mnvMC, dataMCScale, legPos, false, mnvBGMC, mnvBGD, false, true);
 
            if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPos);
            if( WRITE_POT ) mnvPlotter.WriteNorm("POT-Normalized", potX, potY, .04, datapot );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "DIS Signal (BG Subtracted) - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cDIS );
        }

        {
            
            TString cName = Form( "DIS_SignalPlot_dataMCRatio_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cR( cName, cName, 1280, 800 );      
            cR.SetGrayscale( GRAYSCALE );
            mnvMC->GetYaxis()->SetTitle( "Data/MC" );
            TString Yaxis = "Data/MC"; 
             std::string x_label = var.c_str();
	    //mnvPlotter.DrawDataMCRatio( mnvD, mnvMC, dataMCScale, false, true, true, .6, 1.7 );
	    mnvPlotter.DrawDataMCRatio( mnvD, mnvMC, dataMCScale, true, true, 0.0, 2.0, x_label.c_str(), Yaxis );
            mnvPlotter.AddChi2Label(mnvD, mnvMC, dataMCScale, "TL", .06, 0.0, true );
	    if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPosR);
            if( WRITE_POT ) 
                mnvPlotter.WriteNorm("POT-Normalized", "TR", .04, datapot );
            else
                mnvPlotter.WriteNorm("Shape Only", "TR", .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Signal (BG Subtracted)- %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cR );
            
        }

        cout<<"332"<<endl;
        //====== Draw Unfolded Data-MC Overlay Signal ======//
        {
            
            TString cName = Form("DIS_SignalPlot_Unfolded_dataMC_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cDIS( cName, cName, 1280, 800 );
            cDIS.SetGrayscale( GRAYSCALE );
            
            mnvMC_uf->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvMC_uf->GetYaxis()->SetTitle( y_title.c_str() );
            mnvMC_uf->GetXaxis()->SetRangeUser(minValDIS, maxValDIS);
            cout<<"280"<<endl;
            mnvMC_uf->SetMaximum(2000.0);
            mnvD_uf->SetMaximum(2000.0);
            mnvPlotter.DrawDataMCWithErrorBand( mnvD_uf, mnvMC_uf, dataMCScale, legPos, false, NULL, NULL, false, true);
            
            if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPos);
            if( WRITE_POT ) mnvPlotter.WriteNorm("POT-Normalized", potX, potY, .04, datapot );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "DIS Unfolded - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cDIS );
        }
        {
            
            TString cName = Form( "DIS_SignalPlot_Unfolded_dataMCRatio_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cR( cName, cName, 1280, 800 );
            cR.SetGrayscale( GRAYSCALE );
            //mnvMC->GetYaxis()->SetTitle( "Data/MC" );
            TString Yaxis = "Data/MC"; 
             std::string x_label = var.c_str();
            //mnvPlotter.DrawDataMCRatio( mnvD_uf, mnvMC_uf, dataMCScale, false, true, true, .6, 1.7 );
            mnvPlotter.DrawDataMCRatio( mnvD_uf, mnvMC_uf, dataMCScale, false, true, 0.0, 2.0, x_label.c_str(), Yaxis );
            mnvPlotter.AddChi2Label(mnvD_uf, mnvMC_uf, dataMCScale, "TL", .06 , 0.0, true );
            if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPosR);
            if( WRITE_POT )
                mnvPlotter.WriteNorm("POT-Normalized", "TR", .04, datapot );
            else
                mnvPlotter.WriteNorm("Shape Only", "TR", .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Unfolded - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cR );
            
        }



        {
            
            TString cName = Form("DIS_SignalPlot_Unfolded_EffCorr_dataMC_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cDIS( cName, cName, 1280, 800 );
            cDIS.SetGrayscale( GRAYSCALE );
            
            mnvMC_unf_effcorr->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvMC_unf_effcorr->GetYaxis()->SetTitle( y_title.c_str() );
            mnvMC_unf_effcorr->GetXaxis()->SetRangeUser(minValDIS, maxValDIS);
            cout<<"280"<<endl;
            mnvMC_unf_effcorr->SetMaximum(2000.0);
            mnvD_unf_effcorr->SetMaximum(2000.0);
            
            mnvPlotter.DrawDataMCWithErrorBand(mnvD_unf_effcorr, mnvMC_unf_effcorr, dataMCScale, legPos, false, NULL, NULL, false, true);
            
            if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPos);
            if( WRITE_POT ) mnvPlotter.WriteNorm("POT-Normalized", potX, potY, .04, datapot );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "DIS Efficiency Corrected - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cDIS );
        }
        {
            
            TString cName = Form( "DIS_SignalPlot_Unfolded_Eff_Corr_dataMCRatio_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cR( cName, cName, 1280, 800 );
            cR.SetGrayscale( GRAYSCALE );
            TString Yaxis = "Data/MC"; 
             std::string x_label = var.c_str();
            //mnvMC->GetYaxis()->SetTitle( "Data/MC" );
            //mnvPlotter.DrawDataMCRatio( mnvD_unf_effcorr, mnvMC_unf_effcorr, dataMCScale, false, true, true, .6, 1.7 );
            mnvPlotter.DrawDataMCRatio( mnvD_unf_effcorr, mnvMC_unf_effcorr, dataMCScale, false, true, 0.0, 2.0, x_label.c_str(), Yaxis );
            mnvPlotter.AddChi2Label(mnvD_unf_effcorr, mnvMC_unf_effcorr, dataMCScale, "TL", .06, 0.0, true );
            if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPosR);
            if( WRITE_POT )
                mnvPlotter.WriteNorm("POT-Normalized", "TR", .04, datapot );
            else
                mnvPlotter.WriteNorm("Shape Only", "TR", .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Efficiency Corrected - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cR );
            
        }
        
        {
            
            TString cName = Form("DIS_SignalPlot_dataMC_Unfolded_Eff_Corr_Nucleons_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cDIS( cName, cName, 1280, 800 );
            cDIS.SetGrayscale( GRAYSCALE );
            
            mnvMC_unf_effcorr_nucleons->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvMC_unf_effcorr_nucleons->GetYaxis()->SetTitle( y_title.c_str() );
            mnvMC_unf_effcorr_nucleons->GetXaxis()->SetRangeUser(minValDIS, maxValDIS);
            cout<<"280"<<endl;
            mnvMC_unf_effcorr_nucleons->SetMaximum(2000.0);
            mnvD_unf_effcorr_nucleons->SetMaximum(2000.0);
            
            mnvPlotter.DrawDataMCWithErrorBand( mnvD_unf_effcorr_nucleons, mnvMC_unf_effcorr_nucleons, dataMCScale, legPos, false, NULL, NULL, false, true);
            
            if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPos);
            if( WRITE_POT ) mnvPlotter.WriteNorm("POT-Normalized", potX, potY, .04, datapot );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Efficiency Corrected per Nucleon - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cDIS );
        }
        {
            
            TString cName = Form( "DIS_SignalPlot_dataMCRatio_Unfolded_Eff_Corr_Nucleons_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cR( cName, cName, 1280, 800 );
            cR.SetGrayscale( GRAYSCALE );
            TString Yaxis = "Data/MC"; 
             std::string x_label = var.c_str();
            //mnvMC->GetYaxis()->SetTitle( "Data/MC" );
	    //mnvPlotter.DrawDataMCRatio( mnvD_unf_effcorr_nucleons, mnvMC_unf_effcorr_nucleons, dataMCScale, false, true, true, .6, 1.7 );
            mnvPlotter.DrawDataMCRatio( mnvD_unf_effcorr_nucleons, mnvMC_unf_effcorr_nucleons, dataMCScale, false, true, 0.0, 2.0, x_label.c_str(), Yaxis );
            mnvPlotter.AddChi2Label(mnvD_unf_effcorr_nucleons, mnvMC_unf_effcorr_nucleons, dataMCScale, "TL", .06, 0.0, true );
            if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPosR);
            if( WRITE_POT )
                mnvPlotter.WriteNorm("POT-Normalized", "TR", .04, datapot );
            else
                mnvPlotter.WriteNorm("Shape Only", "TR", .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Efficiency Corrected per Nucleon - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cR );
            
        }

        //===============================
        // ERRORS
        //===============================
        mnvPlotter.height_nspaces_per_hist = 1.1;
        mnvPlotter.width_xspace_per_letter = 0.525;
        mnvPlotter.legend_n_columns = 2;
        
  /////      mnvPlotter.axis_maximum  = 0.15;
        string prelimLoc = "TR";
        string legLoc = "TL";
        double prelimX = 0.7; 
        double prelimY = 0.7;
        
  /*      
        {
            
//            mnvPlotter.axis_maximum = err_axis_max;
            TString cName = Form( "DIS_SignalPlot_Data_errSummary_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cErrSummaryDIS( cName, cName, 1280,920 );
            cErrSummaryDIS.SetGrayscale( GRAYSCALE );
            //mnvD->GetXaxis()->SetTitle( x_titleReco.c_str() );
            //mnvD->GetXaxis()->SetTitle( var.c_str() );
            mnvD->GetYaxis()->SetTitle( "Fractional Uncertainty" );
            
            std::string x_label = var.c_str(); 
            PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
            plotter.DrawErrorSummary( mnvD, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
            //plotter.DrawErrorSummary(mnvD, legLoc, true);
            if( 0 )
                mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Error Sources on DIS Signal in %s", targetString.c_str() ), titleSize);
            mnvPlotter.MultiPrint( &cErrSummaryDIS );
          
            //NukeUtils::Get().PrintUncertaintyTables(mnvD, Form("DIS_Signal_Data_%d_%d_%s.tex",  targetID,targetZ,var.c_str()));
        }
        
        {
            
            //            mnvPlotter.axis_maximum = err_axis_max;
            TString cName = Form( "DIS_SignalPlot_Data_errSummary_DetRes_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cErrSummaryDIS( cName, cName, 1280,920 );
            cErrSummaryDIS.SetGrayscale( GRAYSCALE );
            mnvD->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvD->GetYaxis()->SetTitle( "Fractional Uncertainty" );
            
            PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
            plotter.DrawErrorSummary(mnvD, legLoc, true, true, 0.00001 , false , "Detector Res.");
            if( 0 )
                mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Error Sources on DIS Signal in %s", targetString.c_str() ), titleSize);
            mnvPlotter.MultiPrint( &cErrSummaryDIS );
            
            //NukeUtils::Get().PrintUncertaintyTables(mnvD, Form("DIS_Signal_Data_%d_%d_%s.tex",  targetID,targetZ,var.c_str()));
        }
        

        {
            
//            mnvPlotter.axis_maximum = err_axis_max;
            TString cName = Form( "DIS_SignalPlot_MC_errSummary_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cErrSummaryDIS( cName, cName, 1280,920 );
            cErrSummaryDIS.SetGrayscale( GRAYSCALE );
            std::string x_label = var.c_str(); 
            //mnvMC->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvMC->GetYaxis()->SetTitle( "Fractional Uncertainty" );
            
            PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
            //PlotUtils::MnvPlotter plotter(kNukeCCStyle);
            plotter.DrawErrorSummary(mnvMC, legLoc, true );
            plotter.DrawErrorSummary(mnvMC, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
            //mnvPlotter.DrawErrorSummary(mnvMC, legLoc, true );
            //mnvPlotter.DrawErrorSummary(mnvMC, legLoc, true, true, 0.0001, true, " ", " ", var.c_str());
            if( 0 )
                mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Error Sources on DIS Signal in %s", targetString.c_str() ), titleSize);
            mnvPlotter.MultiPrint( &cErrSummaryDIS );
            
             //NukeUtils::Get().PrintUncertaintyTables(mnvMC, Form("DIS_Signal_MC_%d_%d_%s.tex",  targetID,targetZ,var.c_str()));
        }
        
        {
            
//            mnvPlotter.axis_maximum = err_axis_max;
            TString cName = Form( "DIS_SignalPlot_Data_Unfolded_errSummary_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cErrSummaryDIS( cName, cName, 1280,920 );
            cErrSummaryDIS.SetGrayscale( GRAYSCALE );
            //mnvD_uf->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvD_uf->GetYaxis()->SetTitle( "Fractional Uncertainty" );
            std::string x_label = var.c_str(); 
            PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
            plotter.DrawErrorSummary( mnvD_uf, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
            //plotter.DrawErrorSummary(mnvD_uf, legLoc, true );
            if( 0 )
                mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Error Sources on DIS Signal in %s", targetString.c_str() ), titleSize);
            mnvPlotter.MultiPrint( &cErrSummaryDIS );
            
             //NukeUtils::Get().PrintUncertaintyTables(mnvD_uf, Form("DIS_Signal_Data_Unfolded_%d_%d_%s.tex",  targetID,targetZ,var.c_str()));
        }
        
        {
            
            //            mnvPlotter.axis_maximum = err_axis_max;
            TString cName = Form( "DIS_SignalPlot_Data_Unfolded_errSummary_DetRes_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cErrSummaryDIS( cName, cName, 1280,920 );
            cErrSummaryDIS.SetGrayscale( GRAYSCALE );
            mnvD_uf->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvD_uf->GetYaxis()->SetTitle( "Fractional Uncertainty" );
            
            PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
            plotter.DrawErrorSummary(mnvD_uf, legLoc, true, true, 0.00001 , false , "Detector Res.");
            if( 0 )
                mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Error Sources on DIS Signal in %s", targetString.c_str() ), titleSize);
            mnvPlotter.MultiPrint( &cErrSummaryDIS );
            
            //NukeUtils::Get().PrintUncertaintyTables(mnvD_uf, Form("DIS_Signal_Data_Unfolded_%d_%d_%s.tex",  targetID,targetZ,var.c_str()));
        }
        


        {
            
//            mnvPlotter.axis_maximum = err_axis_max;
            TString cName = Form( "DIS_SignalPlot_MC_Unfolded_errSummary_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cErrSummaryDIS( cName, cName, 1280,920 );
            cErrSummaryDIS.SetGrayscale( GRAYSCALE );
            mnvMC_uf->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvMC_uf->GetYaxis()->SetTitle( "Fractional Uncertainty" );
            std::string x_label = var.c_str();  
            PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
            plotter.DrawErrorSummary(mnvMC_uf, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
            //plotter.DrawErrorSummary(mnvMC_uf, legLoc, true );
            if( 0 )
                mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Error Sources on DIS Signal in %s", targetString.c_str() ), titleSize);
            mnvPlotter.MultiPrint( &cErrSummaryDIS );
            
            //NukeUtils::Get().PrintUncertaintyTables(mnvMC_uf, Form("DIS_Signal_MC_Unfolded_%d_%d_%s.tex",  targetID,targetZ,var.c_str()));
        }
       
        {
            
            //            mnvPlotter.axis_maximum = err_axis_max;
            TString cName = Form( "DIS_SignalPlot_MC_Unfolded_errSummary_FSI_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cErrSummaryDIS( cName, cName, 1280,920 );
            cErrSummaryDIS.SetGrayscale( GRAYSCALE );
            mnvMC_uf->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvMC_uf->GetYaxis()->SetTitle( "Fractional Uncertainty" );
            std::string x_label = var.c_str(); 
            PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
            plotter.DrawErrorSummary( mnvMC_uf, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
            //plotter.DrawErrorSummary(mnvMC_uf, legLoc, true, true, 0.00001 , false , "FSI Models");
            if( 0 )
            mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Error Sources on DIS Signal in %s", targetString.c_str() ), titleSize);
            mnvPlotter.MultiPrint( &cErrSummaryDIS );
    
        }
       
        {
            
            //            mnvPlotter.axis_maximum = err_axis_max;
            TString cName = Form( "DIS_SignalPlot_MC_Unfolded_errSummary_DetRes_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cErrSummaryDIS( cName, cName, 1280,920 );
            cErrSummaryDIS.SetGrayscale( GRAYSCALE );
            mnvMC_uf->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvMC_uf->GetYaxis()->SetTitle( "Fractional Uncertainty" );
            
            PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
            plotter.DrawErrorSummary(mnvMC_uf, legLoc, true, true, 0.00001 , false , "Detector Res.");
            if( 0 )
            mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Error Sources on DIS Signal in %s", targetString.c_str() ), titleSize);
            mnvPlotter.MultiPrint( &cErrSummaryDIS );
            
        }
        
        {
            
            //            mnvPlotter.axis_maximum = err_axis_max;
            TString cName = Form( "DIS_SignalPlot_Unfolded_EffCorrected_Data_errSummary_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cErrSummaryDIS( cName, cName, 1280,920 );
            cErrSummaryDIS.SetGrayscale( GRAYSCALE );
            mnvD->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvD->GetYaxis()->SetTitle( "Fractional Uncertainty" );
            std::string x_label = var.c_str(); 
            PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
            //plotter.DrawErrorSummary(mnvD_unf_effcorr, legLoc, true );
            plotter.DrawErrorSummary( mnvD_unf_effcorr, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
            if( 0 )
                mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Error Sources on DIS Signal in %s", targetString.c_str() ), titleSize);
            mnvPlotter.MultiPrint( &cErrSummaryDIS );
            
            //NukeUtils::Get().PrintUncertaintyTables(mnvD, Form("DIS_Signal_Unfolded_EffCorrected_Data_%d_%d_%s.tex",  targetID,targetZ,var.c_str()));
        }
        
        {
            
            //            mnvPlotter.axis_maximum = err_axis_max;
            TString cName = Form( "DIS_SignalPlot_Unfolded_EffCorrected_MC_errSummary_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cErrSummaryDIS( cName, cName, 1280,920 );
            cErrSummaryDIS.SetGrayscale( GRAYSCALE );
            mnvD->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvD->GetYaxis()->SetTitle( "Fractional Uncertainty" );
            
            PlotUtils::MnvPlotter plotter(kCCQENuInclusiveStyle);
            //plotter.DrawErrorSummary(mnvMC_unf_effcorr, legLoc, true );
            plotter.DrawErrorSummary( mnvMC_unf_effcorr, legLoc, true, true, 0.0001, true, "", "", var.c_str()); 
            if( 0 )
                mnvPlotter.WritePreliminary( prelimX, prelimY, .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form("Error Sources on DIS Signal in %s", targetString.c_str() ), titleSize);
            mnvPlotter.MultiPrint( &cErrSummaryDIS );
            
            //NukeUtils::Get().PrintUncertaintyTables(mnvD, Form("DIS_Signal_Unfolded_EffCorrected_MC_%d_%d_%s.tex",  targetID,targetZ,var.c_str()));
        }
*/

        //cleanup
        cout << "   ...cleaning up..." << endl;
	        delete mnvD;
	        delete mnvMC;
        
        mnvPlotter.CleanTmp();



}//variable loop




}

