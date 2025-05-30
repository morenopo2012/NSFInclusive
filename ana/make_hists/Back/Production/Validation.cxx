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
   bool RunCodeWithSystematics = true;
   //bool RunCodeWithSystematics = false;
   TString histFileName1, histFileName2;

   cout<<"******************************************************************************************"<<endl;
   cout<<"                   I am making plots for both before and after tuning!!                   "<<endl;
   cout<<"******************************************************************************************"<<endl;

   if(RunCodeWithSystematics){
       //histFileName1 = "/minerva/data/users/zdar/NukeHists/Mnv210924_git/minervame1L/NukeCCAnne/Hists_PlaneStacks_t3_z26_Nu_Mnv210924_git.root"; 
       histFileName1 = "/minerva/data/users/zdar/NukeHists/Mnv210924_git/minervame1L/MADANNStuff/Hists_PlaneStacks_t2_z26_Nu_Mnv210924_git.root"; 
       histFileName2 = "/minerva/data/users/zdar/NukeHists/Mnv210924_git/minervame1L/NukeCCStuff/Hists_PlaneStacks_t2_z26_Nu_Mnv210924_git.root"; 
       //histFileName2 = "/minerva/data/users/zdar/NukeHists/Mnv210924_git/minervame1L/Hists_PlaneStacks_t3_z26_Nu_Mnv210924_git.root"; 
     }
   else{
       histFileName1 = Form("%s/Hists_PhysicsBackgd_without_SYS_t%d_z%02d_Nu.root", outdir.c_str(), targetID, targetZ);
     } 
   
   cout<<histFileName1<<endl;
   cout<<histFileName2<<endl;

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
   vars.push_back("Ehad");
   vars.push_back("vtxz");
   vars.push_back("planeDNN");
  

   
 
    std::ofstream mcxBjTable;
    
    if(targetID < 10)
        mcxBjTable.open(Form("Event_table_z%02d.tex", targetZ), std::ofstream::out);
    
    else
        mcxBjTable.open(Form("Event_table_z%02d.tex", -1), std::ofstream::out); 
   
  string targetString;  
  if( targetZ == 6 )
    targetString = "Carbon";
  else if( targetZ == 26 )
    targetString = "Iron of target 2";
  if( targetZ == 82 )
    targetString = "Lead of target 4";

   //const string targetString = targetID < 10 ?

   const string suffix = Form( "t%02d_z%02d", targetID, targetZ );
   
   HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(playlist);
   
   TFile *f1 = new TFile( histFileName1,"read" );
   assert(f1);
   mnvPlotter.ApplyAxisStyle( f1 );
 
   TFile *f2 = new TFile( histFileName2,"read" );
   assert(f2);
   mnvPlotter.ApplyAxisStyle( f2 );
 
   
   //TParameter<double> *mcPOT = (TParameter<double>*)f0->Get("MCPOT");
   //TParameter<double> *dataPOT = (TParameter<double>*)f0->Get("DataPOT");
   //double mcpot = mcPOT->GetVal();
   //double datapot = dataPOT->GetVal(); 
   double dataMCScale = 0.483;

   


   std::vector<MnvH1D*> dBGHists, mcBGHists, dHistsBS, mcHistsBS, mcHistsWrongTarget, mcHistsNonDIS;
   std::vector<MnvH1D*> dHists1, dHists2, dHists_unfolded;
   std::vector<MnvH1D*> mcHists, mcHists_unfolded;
   vector< MnvH1D*> hists_Trans_data, hists_Trans_in_Trans, hists_Contin_in_Trans, hists_Signal_in_Trans;
   vector< MnvH1D*> hists_Contin_data, hists_Trans_in_Contin, hists_Contin_in_Contin, hists_Signal_in_Contin;
   vector< MnvH1D*> histos_mc_trans, histos_mc_contin, data_sideband, data_sideband_Contin;
   std::vector<MnvH1D*> effHists;


    const bool useBG = targetID < 10;

     for( vector<string>::iterator ivar = vars.begin(); ivar != vars.end(); ++ivar ){
    
       const string& var = *ivar;
       f1->cd();     
       f2->cd();
     //Save the before background subtracted distributions

       //MnvH1D *d1 = (MnvH1D*)f1->Get(Form("Data_%s_t3_z26", var.c_str()));
       //MnvH1D *d1 = (MnvH1D*)f1->Get(Form("hFe_%sMC_t3_z26", var.c_str()));
       //MnvH1D *d1 = (MnvH1D*)f1->Get(Form("hFe_%sMC_t%02d_z%02d", var.c_str(), targetID, targetZ));
       //MnvH1D *d1 = (MnvH1D*)f1->Get(Form("hPb_%sMC_t4_z82", var.c_str(), targetID, targetZ));
       MnvH1D *d1 = (MnvH1D*)f1->Get(Form("Data_%s_t2_z26", var.c_str(), targetID, targetZ));
       //MnvH1D *d1 = (MnvH1D*)f1->Get(Form("sample_MC_%s_t14_z82", var.c_str()));
        if( !d1 )
            cout << "ERROR: could not get Data1 MnvH1D for variable: " << *ivar << endl;
        dHists1.push_back( d1 );
  
       //MnvH1D *d2 = (MnvH1D*)f2->Get(Form("Data_%s_t3_z26", var.c_str()));
       //MnvH1D *d2 = (MnvH1D*)f2->Get(Form("hFe_%sMC_t3_z26", var.c_str()));
       //MnvH1D *d2 = (MnvH1D*)f2->Get(Form("hFe_%sMC_t%02d_z%02d", var.c_str(), targetID, targetZ));
       MnvH1D *d2 = (MnvH1D*)f2->Get(Form("Data_%s_t2_z26", var.c_str(), targetID, targetZ));
       //MnvH1D *d2 = (MnvH1D*)f2->Get(Form("sample_MC_%s_t14_z82", var.c_str()));
        if( !d2 )
            cout << "ERROR: could not get Data2 MnvH1D for variable: " << *ivar << endl;
        dHists2.push_back( d2 );
    }  


    double sumScaleFactor = 1.0;
    for( unsigned int i = 0; i != vars.size(); ++i )
    {
        const string var = vars[i];
        
  //      cout<<"MC Integral   "<<mcHists[i]->Integral()<<"  "<<mcHists_unfolded[i]->Integral()<<endl;
  //      cout<<"Data Integral "<<dHists[i]->Integral()<<"  "<<dHists_unfolded[i]->Integral()<<endl;
        
        cout << "   Making dataMC plots of " << var << " for (targetID, targetZ, helicity) = "<< targetID << ", " << targetZ << " = " << endl;
        
        
        const double minValDIS = binsDef->GetDISVarMinVal( var );
        const double maxValDIS = binsDef->GetDISVarMaxVal( var );
        
        
        MnvH1D *mnvD  = dHists1[i];
        MnvH1D *mnvMC = dHists2[i];

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

//////////////For before background plots /////////////////////////////



        {
            
            TString cName = Form( "DIS_BeforeBG_dataMCRatio_%s_%s", var.c_str(), suffix.c_str() );
            TCanvas cR( cName, cName, 1280, 800 );      
            cR.SetGrayscale( GRAYSCALE );
            mnvMC->GetYaxis()->SetTitle( "Data/MC" );
            TString Yaxis = "Data MAD / Data NukeCC"; 
             std::string x_label = var.c_str();
	    //mnvPlotter.DrawDataMCRatio( mnvD, mnvMC, dataMCScale, false, true, true, .6, 1.7 );
	    mnvPlotter.DrawDataMCRatio( mnvD, mnvMC, 1.0, true, true, 0.0, 2.0, x_label.c_str(), Yaxis );
            mnvPlotter.AddChi2Label(mnvD, mnvMC, 1.0, "TL", .06, 0.0, true );
	    if( WRITE_PRELIMINARY )
                mnvPlotter.WritePreliminary(prelimPosR);
            //if( WRITE_POT ) 
            //    mnvPlotter.WriteNorm("POT-Normalized", "TR", .04, 1.529 );
            else
                mnvPlotter.WriteNorm("Shape Only", "TR", .04 );
            if(ADD_HISTO_TITLES) mnvPlotter.AddHistoTitle( Form( "Ratio - %s", targetString.c_str() ), titleSize );
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
        


        //cleanup
        cout << "   ...cleaning up..." << endl;
	        delete mnvD;
	        delete mnvMC;
        
        mnvPlotter.CleanTmp();



}//variable loop




}

