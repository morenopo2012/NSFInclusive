#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "include/NukeCCUtilsNSF.h"
#include "include/NukeCC_Binning.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
//#include "syst_common.h"
#include "TArrayD.h"
#include <iostream>
                              
#include "Math/Factory.h"     
#include "Math/Functor.h"     
#include "Math/Minimizer.h"   

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
#include "Cintex/Cintex.h"
#include "../../include/plotting_functions.h"
#ifndef __CINT__
#endif

using namespace NUKECC_ANA;
using namespace PlotUtils;
using namespace std;

const bool do_fits = false;

TH1D* GetVertErrorBandUniverseHist(MnvH1D* mnvh1d, const std::string& errName, int universe); 
void ConstrainPhysBackground(MnvH1D* mnvh1d_trans_data,
                                   MnvH1D* mnvh1d_mc_signal_in_trans,
                                   MnvH1D* mnvh1d_mc_trans_in_trans,
                                   MnvH1D* mnvh1d_mc_contin_in_trans,
                                   MnvH1D* mnvh1d_contin_data,
                                   MnvH1D* mnvh1d_mc_signal_in_contin,
                                   MnvH1D* mnvh1d_mc_trans_in_contin,
                                   MnvH1D* mnvh1d_mc_contin_in_contin,
                                   MnvH1D* mc_trans_in_trans,
                                   MnvH1D* mc_contin_in_trans,
                                   MnvH1D* mc_trans_in_contin,
                                   MnvH1D* mc_contin_in_contin,
                                   MnvH1D* tmp_scale_factor_trans,
                                   MnvH1D* tmp_scale_factor_contin,
                                   double mc_scale, bool onlyCV /*=false*/, bool writeScaleFactor /*=false*/ );

vector<double> CalcPhysScaleFactorMinimizer( TH1D* h1d_data_trans_cv,
                                                  TH1D* h1d_signal_trans_cv,
                                                  TH1D* h1d_trans_trans_cv,
                                                  TH1D* h1d_contin_trans_cv,
                                                  TH1D* h1d_data_contin_cv,
                                                  TH1D* h1d_signal_contin_cv,
                                                  TH1D* h1d_trans_contin_cv,
                                                  TH1D* h1d_contin_contin_cv,
                                                  double mc_scale,
                                                  const std::string& errName);
double getPhysChi2( const double * par);
      TH1D* m_histo_trans_data;
      TH1D* m_histo_sig_trans;
      TH1D* m_histo_trans_trans;
      TH1D* m_histo_contin_trans;
      TH1D* m_histo_contin_data;
      TH1D* m_histo_sig_contin;
      TH1D* m_histo_trans_contin;
      TH1D* m_histo_contin_contin;
 
MnvH1D* TunePhysBackground( MnvH1D *histo, string SB, string outdir, int targetID, int targetZ, string playlist, string var, bool cv_only = false );
MnvH1D* TunePlasticBackground( MnvH1D *tunedHisto, string region, string outdir, int targetZ, string playlist, string var );

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

  TString hname;
  hname= Form("%s/Hists_EventSelection_Physics_WO_sys_t%d_z%02d_Nu_v1_MAD_.root", outdir.c_str(), targetID, targetZ);
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
//  vars.push_back("Ehad");
  vars.push_back("Enu");
  vars.push_back("Emu");

  double datapot = 8.98e+19;//utils->getPOTData(f1);
  double mcpot = 2.00e+20;//utils->getPOTMC(f1); 
  double dataMCScale = datapot/mcpot;
  cout<<"Data POT = "<<datapot<<"   MC POT = "<<mcpot<<"   Data/MC POT = "<<dataMCScale<<endl;

    vector< MnvH1D*> hists_trans_data, hists_contin_data;
    vector< MnvH1D* > hists_true_signal_in_trans,  hists_true_signal_in_contin;
    vector< MnvH1D* > hists_true_trans_in_trans,   hists_true_trans_in_contin;
    vector< MnvH1D* > hists_true_contin_in_contin, hists_true_contin_in_trans;
    vector< MnvH1D* > hists_true_signal_in_signal;
    vector< MnvH1D* > hists_trans_mc_CH_US,  hists_trans_mc_CH_DS;
    vector< MnvH1D* > hists_contin_mc_CH_US, hists_contin_mc_CH_DS;
 
 for( vector<string>::iterator i = vars.begin(); i != vars.end(); ++i )
    {
       const string& var = *i;

       f1->cd();

        hists_true_trans_in_trans.push_back(   (MnvH1D*) f1->Get( Form("hists_trans_in_trans_%s",  var.c_str())));
        hists_true_contin_in_trans.push_back(  (MnvH1D*) f1->Get( Form("hists_contin_in_trans_%s", var.c_str())));
        hists_true_signal_in_trans.push_back(  (MnvH1D*) f1->Get( Form("hists_signal_in_trans_%s", var.c_str())));

        hists_true_trans_in_contin.push_back(  (MnvH1D*) f1->Get( Form("hists_trans_in_contin_%s", var.c_str())));
        hists_true_contin_in_contin.push_back( (MnvH1D*) f1->Get( Form("hists_trans_in_contin_%s", var.c_str())));
        hists_true_signal_in_contin.push_back( (MnvH1D*) f1->Get( Form("hists_trans_in_contin_%s", var.c_str())));

        hists_true_signal_in_signal.push_back( (MnvH1D*) f1->Get( Form("hists_true_signal_%s", var.c_str())));
        hists_trans_mc_CH_US.push_back( (MnvH1D*) f1->Get( Form("hists_trans_mc_CH_US_%s", var.c_str())));
        hists_contin_mc_CH_US.push_back( (MnvH1D*) f1->Get( Form("hists_contin_mc_CH_US_%s", var.c_str())));
        hists_trans_mc_CH_DS.push_back( (MnvH1D*) f1->Get( Form("hists_trans_mc_CH_DS_%s", var.c_str())));
        hists_contin_mc_CH_DS.push_back( (MnvH1D*) f1->Get( Form("hists_contin_mc_CH_DS_%s", var.c_str())));

        hists_trans_data.push_back( (MnvH1D*) f1->Get(Form("selected_data_reco_trans_%s", var.c_str())));
        hists_contin_data.push_back( (MnvH1D*) f1->Get(Form("selected_data_reco_contin_%s", var.c_str())));
cout<<"I am running"<<endl;
 }

   TFile *scaleFactor;
   if(do_fits){
      scaleFactor = new TFile( Form("%s/PhysBGScaleFactors_t%d_z%d_%s_%s.root", outdir.c_str(), targetID, targetZ, playlist.c_str(), getenv("NUKECC_TAG") ), "recreate" );
      assert(scaleFactor);
      }

 for( unsigned int i = 0; i != vars.size(); ++i )
   {
      const string& var = vars[i];
      cout<<"******************************"<<endl;
      cout<<"       variable = "<<var<<endl;
      cout<<"******************************"<<endl;

      string x_title  = "Reconstructed " + binsDef->GetXaxisTitle( var );
      const string y_title  = binsDef->GetYaxisTitle( var );
      string x_title_true = "True " + binsDef->GetXaxisTitle( var );
 
      double x_min = binsDef->GetDISVarMinVal( var );
      double x_max = binsDef->GetDISVarMaxVal( var );
 
      hists_true_trans_in_trans[i]->SetTitle("Transition");
      hists_true_contin_in_trans[i]->SetTitle("Continuum");
      hists_true_signal_in_trans[i]->SetTitle("Signal");

      hists_true_trans_in_contin[i]->SetTitle("Transition");
      hists_true_contin_in_contin[i]->SetTitle("Continuum");
      hists_true_signal_in_contin[i]->SetTitle("Signal");

      MnvH1D *true_signal_in_trans  = (MnvH1D*)hists_true_signal_in_trans[i]->Clone("sig_in_trans");
      MnvH1D *true_signal_in_contin = (MnvH1D*)hists_true_signal_in_contin[i]->Clone("sig_in_contin");

      MnvH1D *true_trans_in_trans_untuned   =  (MnvH1D*)hists_true_trans_in_trans[i]->Clone("trans_in_trans");
      MnvH1D *true_contin_in_trans_untuned  =  (MnvH1D*)hists_true_contin_in_trans[i]->Clone("contin_in_trans");
      MnvH1D *true_trans_in_contin_untuned  =  (MnvH1D*)hists_true_trans_in_contin[i]->Clone("trans_in_contin");
      MnvH1D *true_contin_in_contin_untuned =  (MnvH1D*)hists_true_contin_in_contin[i]->Clone("contin_in_contin");


      MnvH1D *true_trans_in_trans_tuned   = (MnvH1D*)hists_true_trans_in_trans[i]->Clone("trans_in_trans");
      MnvH1D *true_contin_in_trans_tuned  = (MnvH1D*)hists_true_contin_in_trans[i]->Clone("contin_in_trans");
      MnvH1D *true_trans_in_contin_tuned  = (MnvH1D*)hists_true_trans_in_contin[i]->Clone("trans_in_contin");
      MnvH1D *true_contin_in_contin_tuned = (MnvH1D*)hists_true_contin_in_contin[i]->Clone("contin_in_contin");


      TString histName_trans  = "scale_factor_trans_"+ var;
      TString histName_contin = "scale_factor_contin_"+ var;

      int targetZ_tmp, targetID_tmp;
      targetID_tmp = targetID; 
      targetZ_tmp  = targetZ; 

      MnvH1D* tmp_scale_factor_trans   = (MnvH1D*)hists_true_trans_in_trans[i]->Clone(histName_trans);
      MnvH1D* tmp_scale_factor_contin  = (MnvH1D*)hists_true_contin_in_contin[i]->Clone(histName_contin);

            if(targetID < 10){
                vector<MnvH1D*> hists_trans_tot_mc;
                vector<MnvH1D*> hists_contin_tot_mc;
                cout<<"Hello there!   "<<dataMCScale<<endl;
                for( unsigned int i = 0; i != vars.size(); ++i ){
                    //Get the Plastic SB Weighted histograms
                    
                    MnvH1D* tuned_trans_mc_US = (MnvH1D*)hists_trans_mc_CH_US[i]->Clone("tuned_mc_US");
                    MnvH1D* tuned_trans_mc_DS = (MnvH1D*)hists_trans_mc_CH_DS[i]->Clone("tuned_mc_DS");
                    MnvH1D* multiply_trans_mc_US = (MnvH1D*)hists_trans_mc_CH_US[i]->Clone("multiply_mc_US");
                    MnvH1D* multiply_trans_mc_DS = (MnvH1D*)hists_trans_mc_CH_DS[i]->Clone("multiply_mc_DS");
                    
                    multiply_trans_mc_US->Scale(dataMCScale);
                    multiply_trans_mc_DS->Scale(dataMCScale);
                    
                    cout<<"The number of events upstream Transition  "<<multiply_trans_mc_US->Integral()<<endl;
                    cout<<"The number of events downsteam Transition "<<multiply_trans_mc_DS->Integral()<<endl;
                    
                    tuned_trans_mc_US = TunePlasticBackground( multiply_trans_mc_US, "US", outdir.c_str(), targetZ, playlist, vars[i] );
                    tuned_trans_mc_DS = TunePlasticBackground( multiply_trans_mc_DS, "DS", outdir.c_str(), targetZ, playlist, vars[i] );
                    
                    cout<<"The number of events upstream Transition tuned  "<<tuned_trans_mc_US->Integral()<<endl;
                    cout<<"The number of events downsteam Transition tuned "<<tuned_trans_mc_DS->Integral()<<endl;
                    
                    MnvH1D* tuned_contin_mc_US = (MnvH1D*)hists_contin_mc_CH_US[i]->Clone("tuned_mc_US");
                    MnvH1D* tuned_contin_mc_DS = (MnvH1D*)hists_contin_mc_CH_DS[i]->Clone("tuned_mc_DS");
                    MnvH1D* multiply_contin_mc_US = (MnvH1D*)hists_contin_mc_CH_US[i]->Clone("multiply_mc_US");
                    MnvH1D* multiply_contin_mc_DS = (MnvH1D*)hists_contin_mc_CH_DS[i]->Clone("multiply_mc_DS");
                    
                    multiply_contin_mc_US->Scale(dataMCScale);
                    multiply_contin_mc_DS->Scale(dataMCScale);
                    
                    cout<<"The number of events upstream Continuum  "<<multiply_contin_mc_US->Integral()<<endl;
                    cout<<"The number of events downsteam Continuum "<<multiply_contin_mc_DS->Integral()<<endl;
                    
                    tuned_contin_mc_US = TunePlasticBackground( multiply_contin_mc_US, "US", outdir.c_str(), targetZ, playlist, vars[i] );
                    tuned_contin_mc_DS = TunePlasticBackground( multiply_contin_mc_DS, "DS", outdir.c_str(), targetZ, playlist, vars[i] );
                    
                    cout<<"The number of events upstream Continuum tuned  "<<tuned_contin_mc_US->Integral()<<endl;
                    cout<<"The number of events downsteam Continuum tuned "<<tuned_contin_mc_DS->Integral()<<endl;
                    
                    cout<<"I weighted my histograms"<<endl;
                    
                    cout<<"The data integral for the Transition region "<<hists_trans_data[i]->Integral()<<endl;
                    cout<<"The data integral for the Continuum region  "<<hists_contin_data[i]->Integral()<<endl;
                    
                    hists_trans_data[i]->Add(tuned_trans_mc_US,-1.0); //-dataMCScale);
                    hists_trans_data[i]->Add(tuned_trans_mc_DS,-1.0); //-dataMCScale);
                    hists_contin_data[i]->Add(tuned_contin_mc_US,-1.0); //-dataMCScale);
                    hists_contin_data[i]->Add(tuned_contin_mc_DS,-1.0); //-dataMCScale);
                    
                    cout<<"I subtracted my tuned plastic background from my data distributions"<<endl;
                    
                    cout<<"The data integral for the Transition region "<<hists_trans_data[i]->Integral()<<endl;
                    cout<<"The data integral for the Continuum region  "<<hists_contin_data[i]->Integral()<<endl;
                    
                }
            }

   if(do_fits){
      scaleFactor->cd();
      ConstrainPhysBackground(hists_trans_data[i],
                                                true_signal_in_trans,
                                                hists_true_trans_in_trans[i],
                                                hists_true_contin_in_trans[i],
                                                hists_contin_data[i],
                                                true_signal_in_contin,
                                                hists_true_trans_in_contin[i],
                                                hists_true_contin_in_contin[i],
                                                true_trans_in_trans_tuned,
                                                true_contin_in_trans_tuned,
                                                true_trans_in_contin_tuned,
                                                true_contin_in_contin_tuned,
                                                tmp_scale_factor_trans,
                                                tmp_scale_factor_contin,
                                                dataMCScale, false, true);
           
        tmp_scale_factor_trans->Write();
        tmp_scale_factor_contin->Write();
      }
    else{
        true_trans_in_trans_tuned    = TunePhysBackground( true_trans_in_trans_untuned, "trans", outdir.c_str(),  targetID_tmp, targetZ_tmp, playlist, var);
        true_contin_in_trans_tuned   = TunePhysBackground( true_contin_in_trans_untuned, "contin", outdir.c_str(),  targetID_tmp, targetZ_tmp, playlist, var);
        true_trans_in_contin_tuned   = TunePhysBackground( true_trans_in_contin_untuned, "trans", outdir.c_str(),  targetID_tmp, targetZ_tmp, playlist, var);
        true_contin_in_contin_tuned  = TunePhysBackground( true_contin_in_contin_untuned, "contin", outdir.c_str(),  targetID_tmp, targetZ_tmp, playlist, var);

        } 

//===============================ADDING PLOTING STUFF=======================(remove later)
  if(!do_fits)
    {
        PlotUtils::HistFolio<PlotUtils::MnvH1D> hlist_trans_untuned_grouped(Form("trans_untuned_grouped_%s_%s", suffix.c_str(), var.c_str()), true_signal_in_trans);
        hlist_trans_untuned_grouped.AddComponentHist("Trans", true_trans_in_trans_untuned);
        hlist_trans_untuned_grouped.AddComponentHist("Contin", true_contin_in_trans_untuned);
        hlist_trans_untuned_grouped.AddComponentHist("Signal", true_signal_in_trans);

        PlotUtils::HistFolio<PlotUtils::MnvH1D> hlist_contin_untuned_grouped(Form("contin_untuned_grouped_%s_%s", suffix.c_str(), var.c_str()), true_signal_in_contin);
        hlist_contin_untuned_grouped.AddComponentHist("Trans", true_trans_in_contin_untuned);
        hlist_contin_untuned_grouped.AddComponentHist("Contin", true_contin_in_contin_untuned);
        hlist_contin_untuned_grouped.AddComponentHist("Signal", true_signal_in_contin);

        PlotUtils::HistFolio<PlotUtils::MnvH1D> hlist_trans_tuned_grouped(Form("trans_tuned_grouped_%s_%s", suffix.c_str(), var.c_str()), true_signal_in_trans);
        hlist_trans_tuned_grouped.AddComponentHist("Trans", true_trans_in_trans_tuned);
        hlist_trans_tuned_grouped.AddComponentHist("Contin", true_contin_in_trans_tuned);
        hlist_trans_tuned_grouped.AddComponentHist("Signal", true_signal_in_trans);

        PlotUtils::HistFolio<PlotUtils::MnvH1D> hlist_contin_tuned_grouped(Form("contin_tuned_grouped_%s_%s", suffix.c_str(), var.c_str()), true_signal_in_contin);
        hlist_contin_tuned_grouped.AddComponentHist("Trans", true_trans_in_contin_tuned);
        hlist_contin_tuned_grouped.AddComponentHist("Contin", true_contin_in_contin_tuned);
        hlist_contin_tuned_grouped.AddComponentHist("Signal", true_signal_in_contin);

        PlotStacked(hists_trans_data[i], hlist_trans_untuned_grouped.GetHistArray(), dataMCScale, hlist_trans_untuned_grouped.GetName(), hlist_trans_untuned_grouped.GetName());
        PlotStacked(hists_contin_data[i], hlist_contin_untuned_grouped.GetHistArray(), dataMCScale, hlist_contin_untuned_grouped.GetName(), hlist_contin_untuned_grouped.GetName());
        PlotStacked(hists_trans_data[i], hlist_trans_tuned_grouped.GetHistArray(), dataMCScale, hlist_trans_tuned_grouped.GetName(), hlist_trans_tuned_grouped.GetName());
        PlotStacked(hists_contin_data[i], hlist_contin_tuned_grouped.GetHistArray(), dataMCScale, hlist_contin_tuned_grouped.GetName(), hlist_contin_tuned_grouped.GetName());
     }

  TFile* tunedHists;
    tunedHists =  new TFile( Form("%s/TunedPhysicsSidebands_t%d_z%d_%s_%s.root", "."/*outdir.c_str()*/, targetID, targetZ, playlist.c_str(), getenv("NUKECC_TAG") ), "recreate" );

  tunedHists->cd();
   true_trans_in_trans_tuned->Write();
   true_contin_in_trans_tuned->Write();
   true_trans_in_contin_tuned->Write();
   true_contin_in_contin_tuned->Write();

  tunedHists->Close();


     cout<<"done!"<<endl;
  }
 cout<<"done!"<<endl;
 return 0;

}


MnvH1D* TunePhysBackground( MnvH1D *histo, string SB, string outdir, int targetID, int targetZ, string playlist, string var, bool cv_only ){
    cout<<"I'm tuning the Physics background now "<<endl;
    int temp_tgtZ = 0;
    if(targetID > 10)
        temp_tgtZ = 99;
    if(targetID< 10)
        temp_tgtZ = targetZ;
    
    TFile *scaleFactor = new TFile( Form("%s/PhysBGScaleFactors_z%d_%s_%s.root", outdir.c_str(), temp_tgtZ, playlist.c_str(), getenv("NUKECC_TAG") ), "read" );
    cout<<"I'm going to use the following file "<<Form("%s/PhysBGScaleFactor_z%d_%s_%s.root", outdir.c_str(), temp_tgtZ, playlist.c_str(), getenv("NUKECC_TAG") )<<endl;
    assert(scaleFactor);
    cout<<"Got the scale factor file"<<endl;
    scaleFactor->ls();
    
    MnvH1D* h_scaleFactors = (MnvH1D*)scaleFactor->Get( Form("scale_factor_%s_%s", SB.c_str(), var.c_str() ) ); 
    cout<<"The scale factor histogram exists"<<endl;
    cout << "name of scale factors = " << scaleFactor->GetName() << endl;
    int nonzerobin = h_scaleFactors->FindFirstBinAbove(0.);
    MnvH1D* tunedHisto = (MnvH1D*)histo->Clone("tunedHisto");
    MnvH1D* scaleHisto = (MnvH1D*)histo->Clone("scaleHisto");
    scaleHisto->Reset();
    cout<<" I cloned the histogram "<<endl;
    double cv_val = h_scaleFactors->GetBinContent(nonzerobin);
    cout<<"I've got the scale factor stuff for SB "<<SB.c_str()<<"  the CV value is "<<cv_val<<endl;
    
    for(int i=0;i<scaleHisto->GetNbinsX()+2;i++) scaleHisto->SetBinContent(i,cv_val);
   
    tunedHisto->Multiply(tunedHisto,scaleHisto);
    
    return tunedHisto;
  
}

// ===========================================================
// you're entering background constraint zone 
// ============================================================

void ConstrainPhysBackground(MnvH1D* mnvh1d_trans_data,
                                        MnvH1D* mnvh1d_mc_signal_in_trans,
                                        MnvH1D* mnvh1d_mc_trans_in_trans,
                                        MnvH1D* mnvh1d_mc_contin_in_trans,
                                        MnvH1D* mnvh1d_contin_data,
                                        MnvH1D* mnvh1d_mc_signal_in_contin,
                                        MnvH1D* mnvh1d_mc_trans_in_contin,
                                        MnvH1D* mnvh1d_mc_contin_in_contin,
                                        MnvH1D* mc_trans_in_trans,
                                        MnvH1D* mc_contin_in_trans,
                                        MnvH1D* mc_trans_in_contin,
                                        MnvH1D* mc_contin_in_contin,
                                        MnvH1D* tmp_scale_factor_trans,
                                        MnvH1D* tmp_scale_factor_contin,
                                        double mc_scale, bool onlyCV /*=false*/, bool writeScaleFactor /*=false*/ )
{
    std::cout<<"I'm going to constrain the physics background"<<std::endl;
    
    mnvh1d_mc_signal_in_trans->Scale(mc_scale);
    mnvh1d_mc_trans_in_trans->Scale(mc_scale);
    mnvh1d_mc_contin_in_trans->Scale(mc_scale);
    mnvh1d_mc_signal_in_contin->Scale(mc_scale);
    mnvh1d_mc_trans_in_contin->Scale(mc_scale);
    mnvh1d_mc_contin_in_contin->Scale(mc_scale);
    
    mc_trans_in_trans->Scale(mc_scale);
    mc_contin_in_trans->Scale(mc_scale);
    mc_trans_in_contin->Scale(mc_scale);
    mc_contin_in_contin->Scale(mc_scale);
    
    //Get the CV histo
    TH1D* h1d_data_trans_cv = new TH1D(mnvh1d_trans_data->GetCVHistoWithStatError());
    TH1D* h1d_signal_trans_cv = new TH1D(mnvh1d_mc_signal_in_trans->GetCVHistoWithStatError());
    TH1D* h1d_trans_trans_cv = new TH1D(mnvh1d_mc_trans_in_trans->GetCVHistoWithStatError());
    TH1D* h1d_contin_trans_cv = new TH1D(mnvh1d_mc_contin_in_trans->GetCVHistoWithStatError());
    TH1D* h1d_data_contin_cv = new TH1D(mnvh1d_contin_data->GetCVHistoWithStatError());
    TH1D* h1d_signal_contin_cv = new TH1D(mnvh1d_mc_signal_in_contin->GetCVHistoWithStatError());
    TH1D* h1d_trans_contin_cv = new TH1D(mnvh1d_mc_trans_in_contin->GetCVHistoWithStatError());
    TH1D* h1d_contin_contin_cv = new TH1D(mnvh1d_mc_contin_in_contin->GetCVHistoWithStatError());
    
    
    
    std::cout << "Total Trans data:  " << h1d_data_trans_cv->Integral() << std::endl;
    std::cout << "Total Contin data: " << h1d_data_contin_cv->Integral() << std::endl;
    
    vector<double> scales_cv, results_cv;
    // fit tells how much the MC should be scaled to match data
    results_cv = CalcPhysScaleFactorMinimizer(h1d_data_trans_cv,
                                              h1d_signal_trans_cv,
                                              h1d_trans_trans_cv,
                                              h1d_contin_trans_cv,
                                              h1d_data_contin_cv,
                                              h1d_signal_contin_cv,
                                              h1d_trans_contin_cv,
                                              h1d_contin_contin_cv,
                                              mc_scale,
                                              "CV");
    
    //scale only the cv here, universe histograms below
    //scale only the floated histogram;
   
    ((TH1D*)mc_trans_in_trans)->Scale(results_cv[0]);
    ((TH1D*) mc_trans_in_contin)->Scale(results_cv[0]);
    ((TH1D*)tmp_scale_factor_trans)->Divide( ((TH1D*)tmp_scale_factor_trans), ((TH1D*)tmp_scale_factor_trans), results_cv[0], 1.0);
    
   
    ((TH1D*)mc_contin_in_trans)->Scale(results_cv[1]);
    ((TH1D*)mc_contin_in_contin)->Scale(results_cv[1]);
    ((TH1D*)tmp_scale_factor_contin)->Divide( ((TH1D*)tmp_scale_factor_contin), ((TH1D*)tmp_scale_factor_contin), results_cv[1], 1.0);
    
    if( onlyCV ) return;
 
    // vertical error band
    std::vector<std::string> vertNames = mc_trans_in_trans->GetVertErrorBandNames();
    for (std::vector<std::string>::iterator name = vertNames.begin();
         name != vertNames.end(); ++name) {
        std::string errName = *name;
        const unsigned int n_universe = mc_trans_in_trans->GetVertErrorBand(errName)->GetNHists();
        
        //Tune CV in each universe
        if( mc_trans_in_trans->HasVertErrorBand(errName)) {
            MnvVertErrorBand *vertError = mc_trans_in_trans->GetVertErrorBand(errName);
            vertError->Scale(results_cv[0],"",false);
            MnvVertErrorBand *scaleVertError = tmp_scale_factor_trans->GetVertErrorBand(errName);
            scaleVertError->TH1D::Divide( ((TH1D*)tmp_scale_factor_trans), ((TH1D*)tmp_scale_factor_trans), results_cv[0], 1.0);
        }
        if( mc_contin_in_trans->HasVertErrorBand(errName)) {
            MnvVertErrorBand *vertError = mc_contin_in_trans->GetVertErrorBand(errName);
            vertError->Scale(results_cv[1],"",false);
            MnvVertErrorBand *scaleVertError = tmp_scale_factor_contin->GetVertErrorBand(errName);
            scaleVertError->TH1D::Divide( ((TH1D*)tmp_scale_factor_contin), ((TH1D*)tmp_scale_factor_contin), results_cv[1], 1.0);
        }
        if( mc_trans_in_contin->HasVertErrorBand(errName)) {
            MnvVertErrorBand *vertError = mc_trans_in_contin->GetVertErrorBand(errName);
            vertError->Scale(results_cv[0],"",false);
            MnvVertErrorBand *scaleVertError = tmp_scale_factor_trans->GetVertErrorBand(errName);
            scaleVertError->TH1D::Divide( ((TH1D*)tmp_scale_factor_trans), ((TH1D*)tmp_scale_factor_trans), results_cv[0], 1.0);
        }
        if( mc_contin_in_contin->HasVertErrorBand(errName)) {
            MnvVertErrorBand *vertError = mc_contin_in_contin->GetVertErrorBand(errName);
            vertError->Scale(results_cv[1],"",false);
            MnvVertErrorBand *scaleVertError = tmp_scale_factor_contin->GetVertErrorBand(errName);
            scaleVertError->TH1D::Divide( ((TH1D*)tmp_scale_factor_contin), ((TH1D*)tmp_scale_factor_contin), results_cv[1], 1.0);
        }
      
        for (unsigned int i = 0; i < n_universe; ++i) {
            
            std::cout << "\t Uncertainty band name:   " << errName << " universe " << i << std::endl;
            std::cout << "\t\tTotal data:      " << h1d_data_trans_cv->Integral() << std::endl;
            
            TH1D* h1d_signal_trans_universe = GetVertErrorBandUniverseHist(mnvh1d_mc_signal_in_trans, errName, i);
            TH1D* h1d_trans_trans_universe  = GetVertErrorBandUniverseHist(mc_trans_in_trans, errName, i);
            TH1D* h1d_contin_trans_universe = GetVertErrorBandUniverseHist(mc_contin_in_trans, errName, i);
            
            TH1D* h1d_signal_contin_universe = GetVertErrorBandUniverseHist(mnvh1d_mc_signal_in_contin, errName, i);
            TH1D* h1d_trans_contin_universe  = GetVertErrorBandUniverseHist(mc_trans_in_contin, errName, i);
            TH1D* h1d_contin_contin_universe = GetVertErrorBandUniverseHist(mc_contin_in_contin, errName, i);
            
            vector<double> results_universe = CalcPhysScaleFactorMinimizer(h1d_data_trans_cv,
                                                                       h1d_signal_trans_universe,
                                                                       h1d_trans_trans_universe,
                                                                       h1d_contin_trans_universe,
                                                                       h1d_data_contin_cv,
                                                                       h1d_signal_contin_universe,
                                                                       h1d_trans_contin_universe,
                                                                       h1d_contin_contin_universe,
                                                                       mc_scale,
                                                                       errName);
            
            
            assert(results_universe[0] > 0.0);
           
          
            if( mnvh1d_mc_signal_in_contin->HasVertErrorBand(errName)) {
                
                TH1D* universe_trans_trans       = GetVertErrorBandUniverseHist(mc_trans_in_trans, errName, i);
                TH1D* universe_scale_histo_trans = GetVertErrorBandUniverseHist(tmp_scale_factor_trans,   errName, i);
                TH1D* universe_scalefactor_trans = GetVertErrorBandUniverseHist(tmp_scale_factor_trans,   errName, i);
                universe_trans_trans->Scale(results_universe[0]);
                universe_scalefactor_trans->Divide(universe_scale_histo_trans,universe_scale_histo_trans,results_universe[0],1.0);
                
               
                TH1D* universe_contin_contin = GetVertErrorBandUniverseHist(mc_contin_in_contin, errName, i);
                TH1D* universe_scale_histo_contin = GetVertErrorBandUniverseHist(tmp_scale_factor_contin,   errName, i);
                TH1D* universe_scalefactor_contin = GetVertErrorBandUniverseHist(tmp_scale_factor_contin,   errName, i);
                universe_contin_contin->Scale(results_universe[1]);
                universe_scalefactor_contin->Divide(universe_scale_histo_contin,universe_scale_histo_contin,results_universe[1],1.0);
                
                TH1D* universe_trans_contin       = GetVertErrorBandUniverseHist(mc_trans_in_contin, errName, i);
                universe_trans_contin->Scale(results_universe[0]);
               
                TH1D* universe_contin_trans = GetVertErrorBandUniverseHist(mc_contin_in_trans, errName, i);
                universe_contin_trans->Scale(results_universe[1]);
                
                std::cout.precision(8);
                std::cout << "\t\t " << errName << " " << i << " check scale US: " << results_universe[0] << std::endl;
                
            }// end of upstream
            
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


vector<double> CalcPhysScaleFactorMinimizer( TH1D* h1d_data_trans_cv,
                                                       TH1D* h1d_signal_trans_cv,
                                                       TH1D* h1d_trans_trans_cv,
                                                       TH1D* h1d_contin_trans_cv,
                                                       TH1D* h1d_data_contin_cv,
                                                       TH1D* h1d_signal_contin_cv,
                                                       TH1D* h1d_trans_contin_cv,
                                                       TH1D* h1d_contin_contin_cv,
                                                       double mc_scale,
                                                       const std::string& errName)
{
    
    std::cout<<"I'm inside the physics scale minimizer" <<std::endl;
    
    //    cout << "nominal event before fit ( data, signal, other, plastic US, plastic DS ) : " << std::fixed << std::setprecision(2) << h1d_data->Integral() << " & "
    //    << std::fixed << std::setprecision(2) <<  h1d_mc_signal->Integral() << " & "
    //    << std::fixed << std::setprecision(2) <<  h1d_mc_other->Integral() << " & "
    //    << std::fixed << std::setprecision(2) <<  h1d_mc_backgroundUS->Integral() << " & "
    //    << std::fixed << std::setprecision(2) <<  h1d_mc_backgroundDS->Integral() << "\\\\" << endl;
    
    m_histo_trans_data   = h1d_data_trans_cv;
    m_histo_sig_trans    = h1d_signal_trans_cv;
    m_histo_trans_trans  = h1d_trans_trans_cv;
    m_histo_contin_trans = h1d_contin_trans_cv;
    m_histo_contin_data   = h1d_data_contin_cv;
    m_histo_sig_contin    = h1d_signal_contin_cv;
    m_histo_trans_contin  = h1d_trans_contin_cv;
    m_histo_contin_contin = h1d_contin_contin_cv;
    
    std::cout<<"I've got my histograms now" <<std::endl;
    
    
    // initialize stuff
    double scale_factor_trans  = 1.0;
    double scale_factor_contin = 1.0;
    double minChi2 = -1.0;
    
    vector<double> results;
    results.push_back(scale_factor_trans);
    results.push_back(scale_factor_contin);
    results.push_back(minChi2);
    
    // Pick some numbers to be the tolerance and stuff
    ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2");
    fitter->SetMaxFunctionCalls(1000000);
    fitter->SetMaxIterations(100000);
    fitter->SetTolerance(0.01);
    fitter->SetPrintLevel(1);
    
    std::cout<<"Made a Minimizer" <<std::endl;
    
    
    //BkgFitter::getChi2 has 1 parameter, minimizes chi2 for joint fit
    // initialize scale to 1.0 with step size of 0.001
    //    fitter->SetVariable( 0, "scale", 1.0, 0.001 );
   
//    ROOT::Math::Functor lf( this, &getPhysChi2, 2 );
      ROOT::Math::Functor lf( & getPhysChi2, 2 );
    ROOT::Math::Functor functor( lf, 2 );
    fitter->SetFunction( functor );
    fitter->SetLimitedVariable( 0, "trans", 3.0, 0.01, 0.1, 5.0 );
    fitter->SetLimitedVariable( 1, "contin", 3.0, 0.01, 0.1, 5.0 );
    //    fitter->setFitRange(2.0,25.0);

    std::cout<<"Defined my functor " <<std::endl;
    
    // Go!
    fitter->Minimize();
    
    std::cout<<"I minimized!" <<std::endl;
    if( fitter->Status() != 0 ) {
        std::cout << "Womp womp womp... The background fit failed!" << std::endl;
        return results;
    }
    
    const double *bestfit = fitter->X();
    minChi2 = getPhysChi2( bestfit );
    results[0] = bestfit[0]; //scale_factor_trans;
    results[1] = bestfit[1]; //scale_factor_contin;
    results[2] = minChi2;
    
    std::cout<<"The results are "<<results[0]<<"  "<<results[1]<<"  "<<results[2]<<std::endl;
    return results;
    
}


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

MnvH1D*  TunePlasticBackground( MnvH1D *histo, string region, string outdir, int targetZ, string playlist, string var ){
    cout<<"OutDir = "<<outdir<<"  Z = "<<targetZ<<"   playlist = "<<playlist<<"   var = "<<var<<endl;  
    TFile *scaleFactor = new TFile( Form("%s/CHScaleFactors_z%d_%s_%s.root", outdir.c_str(), targetZ, playlist.c_str(), getenv("NUKECC_TAG") ), "read" );
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
    cout<<"done tuning plastic SB"<<endl; 
    return tunedHisto;
    
}
