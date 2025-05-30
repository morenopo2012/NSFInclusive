//include classes
#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "../../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
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
#include "../../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "../../include/Variable.h"
//#include "../../include/Variable_faiza.h"
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
#include "TParameter.h"
#include "../../include/plotting_functions.h"
#ifndef __CINT__
#endif

using namespace NUKECC_ANA;
using namespace PlotUtils;
using namespace std;


void DrawPopulationMigrationHistogram( const TH2D* h_migration, const bool drawAsMatrix, const bool includeFlows, const double binTextSize );
void DrawNormalizedMigrationHistogram( const TH2D* h_migration, const bool drawAsMatrix, const bool includeFlows );

void DrawPopulationMigrationHistogram( const TH2D* h_migration,
                                       const bool drawAsMatrix,
                                       const bool includeFlows,
                                       const double binTextSize )
{
    PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kNukeCCStyle);
    mnvPlotter.SetRedHeatPalette();

    int first_bin = includeFlows ? 0 : 1;
    int last_bin = includeFlows ? h_migration->GetNbinsX()+1 : h_migration->GetNbinsX();
    int nbins = includeFlows ? h_migration->GetNbinsX()+2 : h_migration->GetNbinsX();

    TMatrixD m_migration(nbins, nbins);
    TH2D tmp(*h_migration);
    tmp.Reset();
    for( int y = first_bin; y <= last_bin; ++y ) 
    {
         for( int x = first_bin; x <= last_bin; ++x ) 
         {
              if( includeFlows ) 
                  m_migration[y][x] = h_migration->GetBinContent(x,y); //yeah that's right  y/x
              else
                  m_migration[y-1][x-1] = h_migration->GetBinContent(x,y); //yeah that's right  y/x
         
              tmp.SetBinContent( x, y, h_migration->GetBinContent(x,y));
         }
    }

    if( drawAsMatrix ) 
    {
        tmp = TH2D( m_migration );
        tmp.GetXaxis()->CenterTitle();
        tmp.GetXaxis()->SetTitleOffset( mnvPlotter.axis_title_offset_x );
        tmp.GetYaxis()->SetTitleOffset( mnvPlotter.axis_title_offset_y );
        tmp.GetZaxis()->SetTitleOffset( 0.8 );
        tmp.GetXaxis()->SetNdivisions(100 + h_migration->GetNbinsX() + 2);
        tmp.GetYaxis()->SetNdivisions(100 + h_migration->GetNbinsX() + 2);
        tmp.GetXaxis()->SetTitle( "Reconstructed Bins" );
        tmp.GetYaxis()->SetTitle( "True Bins" );
        tmp.GetZaxis()->SetTitle( "Population in Cell" );
        gStyle->SetPaintTextFormat("2.0f");        
        tmp.SetMarkerSize(binTextSize);
        tmp.DrawCopy("colz text");
    }
    else
    {
        tmp.GetXaxis()->CenterTitle();
        tmp.GetXaxis()->SetTitleOffset( mnvPlotter.axis_title_offset_x );
        tmp.GetYaxis()->SetTitleOffset( mnvPlotter.axis_title_offset_y );
        tmp.GetZaxis()->SetTitleOffset( 0.8 );
        tmp.GetXaxis()->SetNdivisions(509);
        tmp.GetYaxis()->SetNdivisions(509);
        tmp.GetXaxis()->SetTitle( "Reconstructed Bins" );
        tmp.GetYaxis()->SetTitle( "True Bins" );
        tmp.GetZaxis()->SetTitle( "Population in Cell" );
        tmp.DrawCopy("colrz");
    }
}

void DrawNormalizedMigrationHistogram( const TH2D* h_migration,
                                       const bool drawAsMatrix,
                                       const bool includeFlows )
{
    PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kNukeCCStyle);
    mnvPlotter.SetRedHeatPalette();

    int first_bin = includeFlows ? 0 : 1;
    int last_bin = includeFlows ? h_migration->GetNbinsX()+1 : h_migration->GetNbinsX();
    int nbins = includeFlows ? h_migration->GetNbinsX()+2 : h_migration->GetNbinsX();

    TMatrixD m_migration(nbins, nbins);
    TH2D tmp(*h_migration);
    tmp.Reset();
    for( int y = first_bin; y <= last_bin; ++y ) 
    {
         double norm = 0.;
         for( int x = first_bin; x <= last_bin; ++x )
              norm += h_migration->GetBinContent(x,y);

         if( fabs(norm) > 1E-8 ) 
         {
             for( int x = first_bin; x <= last_bin; ++x ) 
             {
                  double percentage = 100 * h_migration->GetBinContent(x,y) / norm;
                  if( includeFlows ) 
                      m_migration[y][x] = percentage; //yeah that's right  y/x
                  else
                      m_migration[y-1][x-1] = percentage; //yeah that's right  y/x
             
                  tmp.SetBinContent( x, y, percentage);
             }
         }
    }

    if( drawAsMatrix ) 
    {
        tmp = TH2D( m_migration );
        tmp.GetXaxis()->CenterTitle();
        tmp.GetXaxis()->SetTitleOffset( mnvPlotter.axis_title_offset_x );
        tmp.GetYaxis()->SetTitleOffset( mnvPlotter.axis_title_offset_y );
        tmp.GetZaxis()->SetTitleOffset( 0.8 );
        tmp.GetXaxis()->SetNdivisions(100 + h_migration->GetNbinsX() + 2);
        tmp.GetYaxis()->SetNdivisions(100 + h_migration->GetNbinsX() + 2);
        tmp.GetZaxis()->SetRangeUser( 1.0, 100. );
        tmp.GetXaxis()->SetTitle( "Reconstructed Bins" );
        tmp.GetYaxis()->SetTitle( "True Bins" );
        tmp.GetZaxis()->SetTitle( "Fraction of Row in Cell" );

        gStyle->SetPaintTextFormat("2.0f");
        tmp.SetMarkerSize(2);
        tmp.DrawCopy("colz text");
    }
    else
    {
        tmp.GetXaxis()->CenterTitle();
        tmp.GetXaxis()->SetTitleOffset( mnvPlotter.axis_title_offset_x );
        tmp.GetYaxis()->SetTitleOffset( mnvPlotter.axis_title_offset_y );
        tmp.GetZaxis()->SetTitleOffset( 0.8 );
        tmp.GetXaxis()->SetNdivisions(509);
        tmp.GetYaxis()->SetNdivisions(509);
        tmp.GetZaxis()->SetRangeUser( 1.0, 100. );
        tmp.GetXaxis()->SetTitle( "Reconstructed Bins" );
        tmp.GetYaxis()->SetTitle( "True Bins" );
        tmp.GetZaxis()->SetTitle( "Fraction of Row in Cell" );
        tmp.DrawCopy("colrz");
    }
}

int main(int argc, char * argv[]){
  //ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if(argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------\
-----"<<std::endl;
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

  TString hname, hname1, hname2;
   if(RunCodeWithSystematics){
  hname2= Form("%s/Hists_Migration_t%d_z%02d_AntiNu_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str());
     }
   else{
  hname2= Form("%s/Hists_Migration_t%d_z%02d_AntiNu_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str() );
     }
  TFile *f1 = new TFile( hname2,"read" );
  //TFile *f2 = new TFile( hname2,"read" );
  cout<<hname1<<endl;

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
  //vars.push_back("Ehad");
  //vars.push_back("x");
  //vars.push_back("y");
  //vars.push_back("Q2");

   HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(playlist);
   //FileType::t_FileType type = utils->GetFileTypeString(std::string fType)
   cout<<"I am here"<<endl;
//-------------------------------

        vector<TH2D*> histsResp;
   
       for( vector<string>::iterator i = vars.begin(); i != vars.end(); ++i )
      {
       const string& var = *i;
       f1->cd();
       
        //TH2D *mnv1  = (TH2D*) f1->Get( Form("selected_Migration_%s",  var.c_str()));
        TH2D *mnv1  = (TH2D*) f1->Get("selected_Migration_Emu_Emu");
        
        cout <<FileType::kAny<< " sample integral = " << mnv1->Integral() << endl;
        if( !mnv1 )
            cout << "ERROR: could not get MC MnvH1D for variable: " << var << endl;
        histsResp.push_back( mnv1 ); 
   
    }
    for( unsigned int ihist = 0; ihist != vars.size(); ++ihist )
    {
         const string& var = vars[ihist];
         const string x_title = binsDef->GetXaxisTitle( var );

         double binTextSize = 1.7;
         if( var == "Emu" || var == "ThetaMu" || var == "Q2" ) binTextSize = 1.5; 
         if( var == "pmu" || var == "pll" || var == "pt" || var == "Enu" ) binTextSize = 1.3; 

         TH2D *hResp = histsResp[ihist];
         if( !hResp )
         {
             cout << "  Can't do anything with a NULL hResp for variable " << var << endl;
             continue;
         }
   
         {
           TString canName = Form( "MigrationPlot_Population_%s_%s_%s", var.c_str(), suffix.c_str(), playlist.c_str() );
           TCanvas c( canName, canName, 1200, 800 );
           c.cd()->SetGridx(false);
           c.cd()->SetGridy(false);
           mnvPlotter.ApplyAxisStyle(hResp);
           DrawPopulationMigrationHistogram(hResp, false, true, binTextSize);
           if( ADD_HISTO_TITLES ) mnvPlotter.AddHistoTitle( Form("%s", x_title.c_str()), .05 );
           mnvPlotter.MultiPrint( &c );
         }

         {
           TString canName = Form( "MigrationPlot_Fraction_%s_%s_%s", var.c_str(), suffix.c_str(), playlist.c_str() );
           TCanvas c( canName, canName, 1200, 800 );
           c.cd()->SetGridx(false);
           c.cd()->SetGridy(false);
           mnvPlotter.ApplyAxisStyle(hResp);
           DrawNormalizedMigrationHistogram(hResp, false, true);
           if( ADD_HISTO_TITLES ) mnvPlotter.AddHistoTitle( Form("%s", x_title.c_str()), .05 );
           mnvPlotter.MultiPrint( &c );
         }

         {
           TString canName = Form( "MigrationPlot_Population_BinByBin_%s_%s_%s", var.c_str(), suffix.c_str(), playlist.c_str() );
           TCanvas c( canName, canName, 1200, 800 );
           c.cd()->SetGridx(false);
           c.cd()->SetGridy(false);
           mnvPlotter.ApplyAxisStyle(hResp);
           DrawPopulationMigrationHistogram(hResp, true, true, binTextSize);
           if( ADD_HISTO_TITLES ) mnvPlotter.AddHistoTitle( Form("%s", x_title.c_str()), .05 );
           mnvPlotter.MultiPrint( &c );
         }
   
         {
           TString canName = Form( "MigrationPlot_Fraction_BinByBin_%s_%s_%s", var.c_str(), suffix.c_str(), playlist.c_str() );
           TCanvas c( canName, canName, 1200, 800 );
           c.cd()->SetGridx(false);
           c.cd()->SetGridy(false);
           mnvPlotter.ApplyAxisStyle(hResp);
           DrawNormalizedMigrationHistogram(hResp, true, true);
           if( ADD_HISTO_TITLES ) mnvPlotter.AddHistoTitle( Form("%s", x_title.c_str()), .05 );
           mnvPlotter.MultiPrint( &c );
         }

    }//loop over variables to plot
   
   
    return 0;
}
