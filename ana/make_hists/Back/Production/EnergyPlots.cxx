#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TArrayD.h" 
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
#include "../../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "../../include/Variable_plasticSB.h"
#include "../../include/Variable_physicsSB.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "../../../NUKECCSRC/ana_common/include/LateralSystematics.h"
#include <iostream>
#include <stdlib.h>
//#include "Cintex/Cintex.h"
#include "../../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"

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
   //bool RunCodeWithSystematics = true;
   bool RunCodeWithSystematics = false;
   TString histFileName, effFileName, eventfile;

   cout<<"******************************************************************************************"<<endl;
   cout<<"                   I am making plots for both before and after tuning!!                   "<<endl;
   cout<<"******************************************************************************************"<<endl;

   if(RunCodeWithSystematics){
       eventfile = Form("%s/Hists_EventSelection_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ ); 
     }
   else{
       eventfile = Form("%s/Hists_EventSelection_t%d_z%02d_AntiNu.root", outdir.c_str(), targetID, targetZ ); 
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
/*   vars.push_back("Emu");
   vars.push_back("Enu");
   vars.push_back("Ehad");
   vars.push_back("x");
   vars.push_back("y");
  */ vars.push_back("planeDNN");
   vars.push_back("vtxz");
  

   
 
//    std::ofstream mcxBjTable;
    
//    if(targetID < 10)
//        mcxBjTable.open(Form("Event_table_z%02d.tex", targetZ), std::ofstream::out);
//    
//    else
//        mcxBjTable.open(Form("Event_table_z%02d.tex", -1), std::ofstream::out); 
   
  string targetString;  
  if( targetZ == 6 )
    targetString = "Carbon";
  else if( targetZ == 26 )
    targetString = "Iron of target 3";
  if( targetZ == 82 )
    targetString = "Lead";

   //const string targetString = targetID < 10 ?

   const string suffix = Form( "t%02d_z%02d", targetID, targetZ );
   
   HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(playlist);
   
   TFile *f1 = new TFile( eventfile,"read" );
   assert(f1);
   mnvPlotter.ApplyAxisStyle( f1 );

   
   TParameter<double> *mcPOT = (TParameter<double>*)f1->Get("MCPOT");
   TParameter<double> *dataPOT = (TParameter<double>*)f1->Get("DataPOT");
   double mcpot = mcPOT->GetVal();
   double datapot = dataPOT->GetVal(); 
   double dataMCScale = datapot/mcpot;
   cout<<"MCPOT = "<<mcpot<<"DataPOT = "<< datapot << "Scale = " << dataMCScale<<endl;

   std::vector<MnvH1D*> dHists;
   std::vector<MnvH1D*> mcHists;


    const bool useBG = targetID < 10;

     for( vector<string>::iterator ivar = vars.begin(); ivar != vars.end(); ++ivar ){
    
       const string& var = *ivar;
       f1->cd();     

       MnvH1D *dsignal = (MnvH1D*)f1->Get(Form("h_data_%s", var.c_str()));
        if( !dsignal )
            cout << "ERROR: could not get Data MnvH1D for variable: " << *ivar << endl;
        dHists.push_back( dsignal );
        
        MnvH1D *mcsignal = (MnvH1D*)f1->Get(Form("h_mc_%s", var.c_str()));
        if( !mcsignal )
            cout << "ERROR: could not get MC MnvH1D for variable: " << *ivar << endl;
        mcHists.push_back( mcsignal );
  

    }  


    double sumScaleFactor = 1.0;
    for( unsigned int i = 0; i != vars.size(); ++i )
    {
        const string var = vars[i];
        
        cout<<"MC Integral   "<<mcHists[i]->Integral()<<endl;
        cout<<"Data Integral "<<dHists[i]->Integral()<<endl;
        
        cout << "   Making dataMC plots of " << var << " for (targetID, targetZ, helicity) = "<< targetID << ", " << targetZ << " = " << endl;
        
        if( ! ( dHists[i] && mcHists[i] ) )
        {
            cout << "  Not enough histograms for variable " << var << endl;
            continue;
        }
        
        

        const double minValDIS = binsDef->GetDISVarMinVal( var );
        const double maxValDIS = binsDef->GetDISVarMaxVal( var );
        

        
        MnvH1D *mnvD  = dHists[i];
        MnvH1D *mnvMC = mcHists[i];
        

        for( int bin=0; bin < mnvMC->GetNbinsX(); bin++ ){
        cout << "bin content, bin error? type? " << bin << ", " << mnvMC->GetBinContent(bin) << ", " << mnvMC->GetBinError(bin) << ", " << FileType::kAny << endl;
        cout << "bin content, bin error? type? " << bin << ", " << mnvD->GetBinContent(bin) << ", " << mnvD->GetBinError(bin) << ", " << FileType::kAny << endl;
        }
        
        
        mnvD->GetXaxis()->SetRangeUser( minValDIS, maxValDIS );
        mnvMC->GetXaxis()->SetRangeUser( minValDIS, maxValDIS );
       
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
        

        //maxY = 2000;
        cout<<"269"<<endl;
        cout<<"MC Integral   "<<mnvMC->Integral()<<endl;
        cout<<"Data Integral "<<mnvD->Integral()<<endl;
        
        

//}
//////////////For before background plots /////////////////////////////

        {
	  TCanvas ctest;
	  MnvH1D* sig2=mnvMC->Clone();
	  cout<<"sig2 nombinwidth "<<sig2->GetNormBinWidth();
	  sig2->Scale(sig2->GetNormBinWidth(), "width");
	  sig2->SaveAs(Form("signal_%s.png", var.c_str() ));
            TString cName = Form("Inclusive_BeforeBG_dataMC_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cDIS( cName, cName, 1280, 800 );
            cDIS.SetGrayscale( GRAYSCALE );
            
            mnvMC->GetXaxis()->SetTitle( x_titleReco.c_str() );
            mnvMC->GetYaxis()->SetTitle( y_title.c_str() );
            mnvMC->GetXaxis()->SetRangeUser(minValDIS, maxValDIS);
            cout<<"280"<<endl;
            mnvMC->SetMaximum(2000.0);
            mnvD->SetMaximum(2000.0);
            

            if(var == "x"){
                mnvD->SetNormBinWidth(0.1);
                mnvMC->SetNormBinWidth(0.1);
            }

            if(var == "Enu"){
                mnvD->SetNormBinWidth(1.0);
                mnvMC->SetNormBinWidth(1.0);
            }

            if(var == "y"){
                mnvD->SetNormBinWidth(0.1);
                mnvMC->SetNormBinWidth(0.1);
            }


	    mnvPlotter.DrawDataMCWithErrorBand( mnvD, mnvMC, dataMCScale, legPos, false, NULL, NULL, false, true);
	    //mnvPlotter.DrawDataMCWithErrorBand( mnvDBS, mnvMCBS, dataMCScale, legPos, false, mnvNonDIS, mnvBGD, false, true);
 
            if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPos);
            if( WRITE_POT ) mnvPlotter.WriteNorm("POT-Normalized", potX, potY, .04, datapot );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Signal - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cDIS );
        }

        {
            
            TString cName = Form( "Inclusive_BeforeBG_dataMCRatio_%s_%s", var.c_str(), suffix.c_str() );
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
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Before BG Subtraction Signal - %s", targetString.c_str() ), titleSize );
            mnvPlotter.MultiPrint( &cR );
            
        }





//////////////////////////////////////////////////////////////////////

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
        

        //cleanup
        cout << "   ...cleaning up..." << endl;
	        delete mnvD;
	        delete mnvMC;
        
        mnvPlotter.CleanTmp();



}//variable loop




}

