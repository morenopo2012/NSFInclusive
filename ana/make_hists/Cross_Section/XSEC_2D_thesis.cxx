//*******Run the code like this: ./CrossSection2DMC .


#include "../../../NUKECCSRC/ana_common/src/NukeCCUtilsNSF.cxx"
#include "../../../NUKECCSRC/ana_common/src/NukeCC_Binning.cxx"

#include "PlotUtils/MnvH2D.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/MinervaUniverse.h"
#include "TFile.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"



using namespace PlotUtils;
using namespace std;

// plot a step in xsec extraction :: PROBLEM - I NEED TO MODIFY THIS FOR PANEL PLOTTING
void Plot(PlotUtils::MnvH1D& hist, const std::string& stepName, const std::string& prefix)
{
  TCanvas can(stepName.c_str());
  hist.GetCVHistoWithError().Clone()->Draw();
  can.Print((prefix + "_" + stepName + ".png").c_str());
 
  // plot uncertainty
  PlotUtils::MnvPlotter plotter;
  plotter.ApplyStyle(PlotUtils::kCCQENuInclusiveStyle);
  plotter.axis_maximum = 0.4;

  plotter.DrawErrorSummary(&hist);
  can.Print((prefix + "_" + stepName + "_uncertaintySummary.png").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Other");
  can.Print((prefix + "_" + stepName + "_otherUncertainties.png").c_str());
}
/*
double GetTotalScatteringCenters(bool isMC)
{
  PlotUtils::TargetUtils trackerInfo;
  double scatteringCenters = trackerInfo.GetTrackerNNucleons( 5980, 8422, isMC, 850 ); //(double minZ, double maxZ, bool isMC, double apothem==850) //PlotUtils::TargetProp::Tracker::Face = 5991.29 //Tracker::Back = 8408.91
  return scatteringCenters;
}
*/

double SetUnfoldingFactor(const std::string& varName, int targetZ) {
double stat_unfolding_factor;

	if(varName=="x_Q2" && targetZ==26)
stat_unfolding_factor=1.46;
else if(varName=="x_Q2" && targetZ==6)
stat_unfolding_factor=1.36;
else if(varName=="x_Q2" && targetZ==82)
stat_unfolding_factor=1.4225;
else if(varName=="W_Q2" && targetZ==82)
stat_unfolding_factor=1.9975;
else if(varName=="W_Q2" && targetZ==26)
stat_unfolding_factor=1.84;
else if(varName=="W_Q2" && targetZ==6)
stat_unfolding_factor=1.76;
else if(varName=="pZmu_pTmu")
stat_unfolding_factor=1e+10;
else if(varName=="x_Q2" && targetZ==99)
stat_unfolding_factor=12;
else if(varName=="W_Q2" && targetZ==99)
stat_unfolding_factor=9.8;
return(stat_unfolding_factor);
}




double GetTotalScatteringCenters(int targetZ, bool isMC)
{
  // TARGET INFO
  PlotUtils::TargetUtils targetInfo;
  double Nucleons;
  
  // Target 1 is generally excluded due to rock muon contamination (in the inclusive analysis)
  if(targetZ == 6){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ); // Target 3
  }
  
  if(targetZ == 26){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 1, targetZ, isMC ) // Target 1
             + targetInfo.GetPassiveTargetNNucleons( 2, targetZ, isMC ) // Target 2
             + targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ); // Target 3
             + targetInfo.GetPassiveTargetNNucleons( 5, targetZ, isMC );// Target 5
  }
  
  if(targetZ == 82){
    	  
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 1, targetZ, isMC )
	     + targetInfo.GetPassiveTargetNNucleons( 2, targetZ, isMC ) // Target 2
             + targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ) // Target 3
             + targetInfo.GetPassiveTargetNNucleons( 4, targetZ, isMC ) // Target 4
             + targetInfo.GetPassiveTargetNNucleons( 5, targetZ, isMC );// Target 5
  }
  if(targetZ > 90 ){
    Nucleons = targetInfo.GetTrackerNNucleons(5980, 8422, isMC, 850);
    //double TargetUtils::GetTrackerNNucleons( double minZ, double maxZ, bool isMC, double apothem /* = 850. */ ) const

  }

  return Nucleons;
}


PlotUtils::MnvH2D* BackgroundSubtract( PlotUtils::MnvH2D* reco, PlotUtils::MnvH2D* bkg, PlotUtils::MnvH2D* data, PlotUtils::MnvH2D* h_bkg_subtracted, double POTNormalization, bool isMC ){
  if(isMC){
    h_bkg_subtracted = (MnvH2D*)reco->Clone("h_background_subtracted_mc");
 // h_bkg_subtracted->Add(bkg,-1);
  }
  else{
   // data->ClearAllErrorBands();
   // data->AddMissingErrorBandsAndFillWithCV(*reco);
    h_bkg_subtracted = (MnvH2D*)data->Clone("h_background_subtracted_data");
   // PlotUtils::MnvH2D* h_scaled = (MnvH2D*)bkg->Clone("h_background_mc_scale");
   // h_scaled->Scale(POTNormalization);
   // h_bkg_subtracted->Add(h_scaled,-1);
  }
  return h_bkg_subtracted;
}




PlotUtils::MnvH2D* UnfoldHist( PlotUtils::MnvH2D* h_data_nobck, PlotUtils::MnvH2D* h_mc_nobck, PlotUtils::MnvH2D* h_reco, PlotUtils::MnvH2D* h_generated, PlotUtils::MnvH2D* h_migration, int num_iter, bool returnDataUnfolded,int targetID, int targetZ,std::string varName)
{
  /*
    TParameter<double> *pot_mc = (TParameter<double>*)f2.Get("MCPOT");
    TParameter<double> *pot_data = (TParameter<double>*)f2.Get("DataPOT");
    double MCPOT = pot_mc->GetVal();
    double DataPOT = pot_data->GetVal();
    double pot_scale = DataPOT/MCPOT;
    h_migration->Scale(pot_scale);
    h_reco->Scale(pot_scale);
    h_generated->Scale(pot_scale);*/
 // POT normalization done later 
 
  MnvH2D *h_data_unfolded = NULL, *h_mc_unfolded = NULL;

  MinervaUnfold::MnvUnfold unfold;
  bool data_unfolded;
  bool mc_unfolded;
  if (returnDataUnfolded)  data_unfolded = unfold.UnfoldHisto2D(h_data_unfolded, h_migration, h_reco, h_generated, h_data_nobck, num_iter, true, false); // last input variable = true makes sure that each universe of the distribution is getting unfolded with its corresponding universe in the migration matrix
  else  mc_unfolded = unfold.UnfoldHisto2D(h_mc_unfolded, h_migration, h_reco, h_generated, h_mc_nobck, num_iter, true, false);

  if(!data_unfolded || !mc_unfolded){
    cout << "Unfolding failed for either data or MC" << endl;
  }

  TCanvas c;
  
  if(data_unfolded){
    cout << "Unfolded data sucessfull for either data or MC" <<endl;
    h_data_unfolded->Draw("COLZ");
    c.SaveAs("data_unfolded.png");
  }

  if(mc_unfolded){
    cout << "Unfolded mc sucessfull for either data or MC" <<endl;
    //h_mc_unfolded->Draw("COLZ");
    //c.SaveAs("mc_unfolded.png");
  }




  // NOW we have to get the covariance of the unfolding which means doing the unfolding again
  // only need this on the data as the MC is used as CV only
  cout << "Getting the covariance of the unfolding" << endl;
  if (returnDataUnfolded){
  TH2D* hUnfoldedDummy=new TH2D(h_data_unfolded->GetCVHistoWithStatError());
  TH2D* hMigrationDummy=new TH2D(h_migration->GetCVHistoWithStatError());
  TH2D* hRecoDummy=new TH2D(h_reco->GetCVHistoWithStatError());
  TH2D* hTruthDummy=new TH2D(h_generated->GetCVHistoWithStatError());
  TH2D* hBGSubDataDummy=new TH2D(h_data_nobck->GetCVHistoWithStatError());

  TMatrixD unfoldingCovMatrixOrig;
  unfold.UnfoldHisto2D(hUnfoldedDummy, unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy, num_iter);

  // There's a bug in RooUnfold that's making it return covariance matriced with two extra bins. Kill them here, wiht a check. Conveniently, this bug was being hidden by an offsetting bug in MnvH2D, which is now fixed (Dan)
  int correctNbins=hUnfoldedDummy->fN;
  int  matrixRows=unfoldingCovMatrixOrig.GetNrows();
  if(correctNbins!=matrixRows){
    cout <<
     "*************************************************************************" << endl;
    cout << "* Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << endl;
   cout <<
     "*************************************************************************" << endl;
   unfoldingCovMatrixOrig.ResizeTo(correctNbins, correctNbins);
  }

  for(int i=0; i<unfoldingCovMatrixOrig.GetNrows(); ++i)
    unfoldingCovMatrixOrig(i,i)=0; //subtract off the diagonal errors on the covariance matrix b/c those are already included elsewhere, we don't want to double count the stat error

  delete hUnfoldedDummy;
  delete hMigrationDummy;
  delete hRecoDummy;
  delete hTruthDummy;
  delete hBGSubDataDummy;

  h_data_unfolded->PushCovMatrix("unfoldingCov", unfoldingCovMatrixOrig);
 double factor = SetUnfoldingFactor(varName, targetZ);
  h_data_unfolded->ModifyStatisticalUnc(factor, "unfoldingCov"); //correction factor determined using additional statistical uncertainties to record uncertainty due to finite MC stat. for 6 ITERATIONS // April 11/2023

  c.Clear();
  unfoldingCovMatrixOrig.Draw("COLZ");
  unfoldingCovMatrixOrig.Print();
  c.SaveAs("unfoldingMatrix.png");

  return h_data_unfolded;
  }
  else return h_mc_unfolded;

} //end UnfoldHisto

// only really have to get the covariance and correlation matrices when someone asks for them
// usually for covariance and correlation matrices want to look at the total error
// Anezka's code: Write covariance and correlation matrices
// Inspired by Gonzalo Diaz's function
void writeCovAndCorrMatrices(const auto&var, PlotUtils::MnvH2D* unfolded_histo, TFile* fout, const auto&prefix)
{
  fout->cd();
 
  // Unfolding covariance matrix
  // It's preferably that, once extracted, you can check unfold_cov_matrix.
  // If your unfolding has been done correctly, it will be a symmetric matrix with no entries in the diagonal. This is because the entries on the diagonal were set to 0 so as not to avoid double counting (since the diagonal values get taken care of elsewhere). 
  // Get syst. error matrix
  // 1st bool: Show as fractional error?
  // 2nd bool: area-norm covariance?
  TMatrixD unfold_cov_matrix = unfolded_histo->GetSysErrorMatrix("unfoldingCov", false, false);
  TH2D* h_unfold_cov_matrix = new TH2D(unfold_cov_matrix); //casting a TMatrix to a TH2D
  h_unfold_cov_matrix->Write(Form("unfold_cov_matrix_%s_%s", prefix.c_str(), var.c_str()));

  // Statistical error matrix
  // bool: show as fractional error?
  TMatrixD stat_err_matrix = unfolded_histo->GetStatErrorMatrix(false);
  TH2D* h_stat_err_matrix = new TH2D(stat_err_matrix);
  h_stat_err_matrix->Write(Form("stat_err_matrix_%s_%s", prefix.c_str(), var.c_str()));

  // Statistical covariance matrix
  // = sum of the unfolding covariance and the stat-only error: not fully diagonal
  // includes the statistical errors of each bin AND the effect of the unfolding
  // between different bins, which in fact adds non-diagonal terms
  TMatrixD stat_cov_matrix = unfold_cov_matrix + stat_err_matrix;
  TH2D* h_stat_cov_matrix = new TH2D(stat_cov_matrix);
  h_stat_cov_matrix->Write(Form("stat_cov_matrix_%s_%s", prefix.c_str(), var.c_str()));

  // Statistical correlation matrix
  const int size = stat_cov_matrix.GetNrows();
  TMatrixD stat_corr_matrix(size, size);
  for ( int x = 0; x < size; ++x) { // loop over columns
    for ( int y = 0; y < size; ++y) { // loop over rows
      stat_corr_matrix[x][y] = (stat_cov_matrix[x][x] == 0.0 || stat_cov_matrix[y][y]) ? 
                                                         0.0 : stat_cov_matrix[x][y] /
                                                         std::sqrt(stat_cov_matrix[x][x]*stat_cov_matrix[y][y]);
    }
  }
  TH2D* h_stat_corr_matrix = new TH2D(stat_corr_matrix);
  h_stat_corr_matrix->Write(Form("stat_corr_matrix_%s_%s", prefix.c_str(), var.c_str()));

  // Total covariance and correlations
  // write total error matrix
  // 1st bool: Include stat. error?
  // 2nd bool: show as fractional error?
  // 3rd bool: area-norm covariance?
  TH2D* h_cov_matrix = new TH2D(unfolded_histo->GetTotalErrorMatrix(true, false, false));
  h_cov_matrix->Write(Form("total_cov_matrix_%s_%s", prefix.c_str(), var.c_str()));

  // Write total correlation matrix
  // 1st bool: Area-norm covariance?
  TH2D* h_corr_matrix = new TH2D(unfolded_histo->GetTotalCorrelationMatrix(false));
  // 2nd bool: Include stat. error?
//  TH2D* h_corr_matrix = new TH2D(unfolded_histo->GetTotalCorrelationMatrixTH2D(false, true));
  
  h_corr_matrix->Write(Form("total_corr_matrix_%s_%s", prefix.c_str(), var.c_str()));
 
  delete h_unfold_cov_matrix;
  delete h_stat_err_matrix;
  delete h_stat_cov_matrix;
  delete h_stat_corr_matrix;
  delete h_cov_matrix;
  delete h_corr_matrix;

} //end of WriteCovAndCorrMatrices func

//The final step of xsec extraction: normalize by flux, bin width, POT, and numeber of targets
// Differential xsec
PlotUtils::MnvH2D* normalize(PlotUtils::MnvH2D* efficiencyCorrected, PlotUtils::MnvH2D* fluxIntegral, const double nNucleons, const double POT)
{
  efficiencyCorrected->Divide(efficiencyCorrected, fluxIntegral); 
  efficiencyCorrected->Scale(1./nNucleons/POT);
  efficiencyCorrected->Scale(1.e4); //Flux histogram is in m^-2 but convention is to report in cm^2
  efficiencyCorrected->Scale(1., "width"); //convention is to store histograms not bin width normalized
  return efficiencyCorrected;  
}

// total xsec -> difference in flux
PlotUtils::MnvH2D* normalizeTotal(PlotUtils::MnvH2D* efficiencyCorrected, PlotUtils::MnvH2D* fluxRebinned, const double nNucleons, const double POT)
{
  efficiencyCorrected->Divide(efficiencyCorrected, fluxRebinned); //PROBLEM: dividing 2D hist by 1D hist
  efficiencyCorrected->Scale(1./nNucleons/POT);
  efficiencyCorrected->Scale(1.e4); //Flux histogram is in m^-2 but convention is to report in cm^2
  return efficiencyCorrected;
}

// ----------------------------------------------------------------------------------------------- //
// ----------------------------------- MAIN ------------------------------------------------------ //
// ----------------------------------------------------------------------------------------------- //
int main(int argc, char *argv[])
{
  TH1::AddDirectory(false); //sets a global switch disabling the referencing so that when a file is closed, all histograms in memory associated with this file are not automatically deleted
  if(argc==1){
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "MACROS HELP:\n\n" <<
    "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n" <<
    "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
//  "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n" <<
//  "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n"
    std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    return 0;
  }
  TString dir(argv[1]);
int targetID = atoi(argv[2]);
int targetZ = atoi(argv[3]);
//  int targetID = 3; int targetZ = 6;
std::string varName=argv[4];
  const std::string plist_string("minervame6A");

  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);

  bool RunCodeWithSystematics = false;
  TString eventLoop, migration, efficiency;
  if (RunCodeWithSystematics){
	  if(targetZ==99)
	  { eventLoop=Form("Bkg_subtracted_tracker_thesis_%s_fullp4_tuneV4_t%d_z%d.root",varName.c_str(),targetID,targetZ);
            migration= Form("Migration_thesis_fullp4_tuneV4_t%d_z%d_sys.root",targetID,targetZ);
	    efficiency=Form("Efficiency_fullp4_tracker_thesis_tuneV4_t%d_z%d_sys.root",targetID,targetZ);
	  }
	  else{
	  eventLoop = Form("Bkg_Sub_thesis_fullp4_%s_tuneV4_SF_pTmu_t%d_z%02d_sys.root",varName.c_str(),targetID,targetZ);
    migration = Form("Migration_fullp4_ns_targetCombined_thesis_tuneV4_t%d_z%02d_sys.root",targetID,targetZ);
    efficiency =Form("Efficiency_fullp4_ns_targetCombined_thesis_tuneV4_t3_z06_sys.root",targetID,targetZ);}
  }
     else{
    eventLoop ="../Hists_EventSelection_t99_z99_Nu_minervame1L.root";
    migration = "../Migration/Migration_minervame1L_t99_z99_nosys.root";
    efficiency ="../Migration/Efficiency_1D_minervame1L_t99_z99_nosys.root";

  }
 
  TFile *fEventLoop = new TFile(eventLoop, "read");
  TFile *fMigration = new TFile(migration, "read");
  TFile *fEfficiency = new TFile(efficiency, "read");


  if (!fEventLoop || !fEventLoop->IsOpen() || fEventLoop->IsZombie()) {
    std::cerr << "Error opening EventLoop file: " << eventLoop << std::endl;
}

// Check Migration file
if (!fMigration || !fMigration->IsOpen() || fMigration->IsZombie()) {
    std::cerr << "Error opening Migration file: " << migration << std::endl;
}

// Check Efficiency file
if (!fEfficiency || !fEfficiency->IsOpen() || fEfficiency->IsZombie()) {
    std::cerr << "Error opening Efficiency file: " << efficiency << std::endl;
}



  NukeCCUtilsNSF *utils = new NukeCCUtilsNSF(plist_string);
  NukeCC_Binning *binsDef = new NukeCC_Binning();
  HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);

  // to iterate over variables
  std::vector<string> vars;
  vars.push_back("Eavailable_q3");
  //vars.push_back("x_Q2");
  //vars.push_back("W_Q2");

  // to iterate over mc and data
  std::vector<string> prefixes;
  prefixes.push_back("mc");
  prefixes.push_back("data");

  // create the output file
  TFile *fXSec = new TFile(Form("XSec_thesis_nuclearflux_%s_t%d_z%02d.root",varName.c_str(),targetID,targetZ), "recreate");


  // read in the POT info
  TParameter<double> *MCPOT = (TParameter<double>*)fEventLoop->Get("MCPOT");
  TParameter<double> *DataPOT = (TParameter<double>*)fEventLoop->Get("DataPOT");
  double mcPOT = MCPOT->GetVal();
  double dataPOT = DataPOT->GetVal();
//    double mcPOT = 4.96197407132E+21;
//    double dataPOT = 1.11707110844e+21;
  // total POT of all RHC
//  double mcPOT = 6.34384E+20;
//  double dataPOT = 1.51022E+20;

  fXSec->cd();
  //write POT to output file
  auto mcPOTOut = new TParameter<double>("MCPOT", mcPOT);
  auto dataPOTOut = new TParameter<double>("DataPOT",dataPOT);
  dataPOTOut->Write();
  mcPOTOut->Write();

  double POTscale = dataPOT/mcPOT;
  cout << "MC: " << mcPOT << "\n" << "Data: " << dataPOT << "\nscale: " << POTscale << endl;

  for (const auto&var : vars){
    //auto migration = dynamic_cast<MnvH2D*>(fMigration->Get(Form("response2d_new_%s_migration", varName.c_str())));
    //auto migration = dynamic_cast<MnvH2D*>(fMigration->Get(Form("selected_mc_response2d_%s_migration", varName.c_str())));
    //auto effNum = dynamic_cast<MnvH2D*>(fEfficiency->Get(Form("h_mc_%s", varName.c_str())));
    //auto effDen = dynamic_cast<MnvH2D*>(fEfficiency->Get(Form("h_truth_%s", varName.c_str())));
    auto migration = dynamic_cast<MnvH2D*>(fMigration->Get("selected_mc_response2d_Eavailable_q3_migration"));
    auto effNum = dynamic_cast<MnvH2D*>(fEfficiency->Get("h_mc_Eavailable_q3"));
    auto effDen = dynamic_cast<MnvH2D*>(fEfficiency->Get("h_truth_Eavailable_q3"));



if (!migration || migration->GetEntries() == 0) {
    std::cerr << "Migration histogram missing or empty!" << std::endl;
}

if (!effNum || effNum->GetEntries() == 0) {
    std::cerr << "Efficiency numerator histogram missing or empty!" << std::endl;
}

if (!effDen || effDen->GetEntries() == 0) {
    std::cerr << "Efficiency denominator histogram missing or empty!" << std::endl;
}


      //effDen->Clone()->Write(Form("Eff_Denominator_%s", varName.c_str()));
      //effNum->Clone()->Write(Form("Eff_Numerator_%s", varName.c_str()));
      effDen->Clone()->Write("Eff_Denominator_Eavailable_q3");
      effNum->Clone()->Write("Eff_Numerator_Eavailable_q3");
    auto simEventRate = effDen->Clone(); //make a copy for later 

 /*   auto reco = dynamic_cast<MnvH2D*>(fEventLoop->Get(Form("h_mc_%s", var.c_str())));
    auto signal = dynamic_cast<MnvH2D*>(fEventLoop->Get(Form("h_mc_signal_%s", var.c_str())));
    MnvH2D* bkg = (MnvH2D*)reco->Clone("bkg");
    bkg->Add(signal, -1);
    bkg->Clone()->Write(Form("total_bkg_%s", var.c_str()));
    auto data = dynamic_cast<MnvH2D*>(fEventLoop->Get(Form("h_data_%s", var.c_str())));*/
  
    //auto reco = dynamic_cast<MnvH2D*>(fEventLoop->Get("h_mc_subtracted"));//added signal instead of h_mc_subtracted 18March,2024
    //auto bkg = dynamic_cast<MnvH2D*>(fEventLoop->Get("h_mc_bkg"));
    //auto data = dynamic_cast<MnvH2D*>(fEventLoop->Get("h_data_subtracted"));
    auto reco = dynamic_cast<MnvH2D*>(fEventLoop->Get("h_mc_Eavailable_q3"));//added signal instead of h_mc_subtracted 18March,2024
    //auto bkg = dynamic_cast<MnvH2D*>(fEventLoop->Get("h_mc_bkg"));
    auto data = dynamic_cast<MnvH2D*>(fEventLoop->Get("h_data_Eavailable_q3"));
    fXSec->cd();

    //PlotUtils::MnvH2D* h_background_subtracted_mc = BackgroundSubtract(reco, bkg, data, h_background_subtracted_mc, POTscale, true);
    //h_background_subtracted_mc->Clone()->Write(Form("Background_subtracted_mc_%s", varName.c_str()));
    //PlotUtils::MnvH2D* h_background_subtracted_data = BackgroundSubtract(reco, bkg, data, h_background_subtracted_data, POTscale, false);
    //h_background_subtracted_data->Clone()->Write(Form("Background_subtracted_data_%s", varName.c_str()));
    
    PlotUtils::MnvH2D* h_background_subtracted_mc = reco;
    h_background_subtracted_mc->Clone()->Write("Background_subtracted_mc_Eavailable_q3");
    PlotUtils::MnvH2D* h_background_subtracted_data = data;
    h_background_subtracted_data->Clone()->Write("Background_subtracted_data_Eavailable_q3");

    //MnvH2D* h_reco=(MnvH2D*)fMigration->Get(Form("%s_reco",varName.c_str())); //branch filled in w/ reco //Prameet commented
    //MnvH2D* h_generated=(MnvH2D*)fMigration->Get(Form("%s_truth",varName.c_str())); //branch filled in w/ true //Prameet commented
    MnvH2D* h_reco=(MnvH2D*)fMigration->Get("h_mc_reco_Eavailable_q3"); //branch filled in w/ reco //Prameet commented
    MnvH2D* h_generated=(MnvH2D*)fMigration->Get("h_mc_true_Eavailable_q3"); //branch filled in w/ true //Prameet commented

    migration->Clone()->Write("selected_mc_response2d_Eavailable_q3_migration");
  
    // added in Aug. 2
    h_reco->Clone()->Write("reco_filled_Eavailable_q3");    //prameet commented 4 materr
    h_generated->Clone()->Write("truth_filled_Eavailable_q3");//prameet commented 4 mar 
    /////// 

 int n_iter;
 if (varName=="x_Q2")
	 n_iter=10;
 else if (varName=="W_Q2" && targetZ!=99)
	 n_iter=5;
 else if(varName=="W_Q2" && targetZ==99)
         n_iter=7;
 else if(varName=="pZmu_pTmu")
	 n_iter=10;
 else if(varName=="Eavailable_q3")
	 n_iter=5;
    // d'Agostini unfolding
    auto unfolded_mc = UnfoldHist(h_background_subtracted_data, h_background_subtracted_mc, h_reco, h_generated, migration,n_iter, false, targetID,targetZ,varName); // passing in 6 iterations here
  
    //if(!unfolded) throw std::runtime_error(std::string("Failed to unfold ") + h_background_subtracted_data->GetName() + " using " + migration->GetName());

    unfolded_mc->Clone()->Write("bkg_subtracted_mc_unfolded_Eavailable_q3");

    auto unfolded = UnfoldHist(h_background_subtracted_data, h_background_subtracted_mc, h_reco, h_generated, migration, n_iter ,true,targetID,targetZ,varName); 
    unfolded->Clone()->Write("bkg_subtracted_data_unfolded_Eavailable_q3");

    std::cout << "Survived writing the unfolded histogram.\n" << std::flush;
    auto templateHist = effNum->Clone(); //needed for passing into GetIntegratedFluxReweighted for copying the binning and systematics
    effNum->Divide(effNum, effDen, 1.0, 1.0, "B"); //handles systematics correctly
    effNum->Clone()->Write("Efficiency_Eavailable_q3");

    // unfolded, efficiency corrected
    unfolded->Divide(unfolded, effNum);
    unfolded->Clone()->Write(Form("Unfolded_efficiencyCorrected_%s", varName.c_str()));
    unfolded_mc->Divide(unfolded_mc, effNum);
    unfolded_mc->Clone()->Write(Form("Unfolded_mc_efficiencyCorrected_%s", varName.c_str()));

    std::cout << "Efficiency corrected!" << std::endl;         

   string material;
      if(targetZ == 6) material = "carbon";
      else if(targetZ == 26) material = "iron";
      else if(targetZ == 82) material = "lead";
      else material = "tracker"; 

  
  // get the integrated flux needed at the end of a differential cross section extraction
  // (Matches binning of input histogram)
  const bool useMuonCorrelations =true;// true;
  int n_flux_universes = PlotUtils::MinervaUniverse::GetNFluxUniverses();
  int nu_pdg = PlotUtils::MinervaUniverse::GetAnalysisNuPDG();
  const bool use_nue_constraint = true;
  double min_energy = 0.0;
  double max_energy = 120.; //Dan uses 100, removing the neutrino energy cut
  const std::string project_dir = "targets_2345_jointNueIMD";
  
  auto& frw = PlotUtils::flux_reweighter(plist_string, nu_pdg, use_nue_constraint, n_flux_universes);
  
//  auto fluxIntegrated = frw.GetIntegratedFluxReweighted(nu_pdg, templateHist, min_energy, max_energy, useMuonCorrelations);
   
   auto fluxIntegrated = frw.GetIntegratedTargetFlux(nu_pdg,material, templateHist, min_energy, max_energy, project_dir); 


  fluxIntegrated->Clone()->Write(Form("FluxIntegrated_%s", varName.c_str())); 
 
  double nNucleons = GetTotalScatteringCenters(targetZ, true); //Tutorial: Always use MC number of nucleons for cross section. Dan: use the same truth fiducial volume for all extractions. The acceptance correction corrects data back to this fiducial even if the reco fiducial cut is different.
fXSec->cd();
  auto crossSection = normalize(unfolded, fluxIntegrated, nNucleons, dataPOT);  
  crossSection->Clone()->Write(Form("CrossSection_%s", varName.c_str()));
fXSec->cd();
  auto crossSection_mc = normalize(unfolded_mc, fluxIntegrated, nNucleons, mcPOT);
  crossSection_mc->Clone()->Write(Form("CrossSectionMC_%s", varName.c_str()));
fXSec->cd();

  //Write a "simulated cross section" to compare to the data I just extracted
  // If this analysis passed its closure test, this should be the same cross section as what GENIEXSecExtract would produce
  normalize(simEventRate, fluxIntegrated, nNucleons, mcPOT);
  simEventRate->Clone()->Write(Form("simulatedCrossSection_%s", varName.c_str()));
  fXSec->cd();
  std::cout << "Data POT: " << dataPOT << std::endl;
  std::cout << "nNucleons: " << nNucleons << std::endl;
  std::cout << "Completed." << std::endl;

fXSec->Close();
delete fXSec;


}
} //end main



