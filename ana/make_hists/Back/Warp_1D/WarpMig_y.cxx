//==============================================================================
// Loop entries, make cuts, fill histograms.
// * Uses the New Systematics Framework and "Universe" objects.
// * loop universes, make cuts and fill histograms with the correct lateral
// shifts and weights for each universe.
// * TChain --> PlotUtils::ChainWrapper.
// * MnvHXD --> PlotUtils::HistWrapper.
// * Genie, flux, non-resonant pion, and some detector systematics calculated.
//==============================================================================

//#include "include/CommonIncludes.h"
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
#include "include/NukeCCUtilsNSF.h"
#include "include/NukeCC_Cuts.h"



#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
//#include "PlotUtils/MnvTuneSystematics.h"



// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we
// shield it.
#ifndef __CINT__
#include "../../include/plotting_functions.h"
#endif
#include "PlotUtils/MacroUtil.h" 
//using namespace globalV;
using namespace NUKECC_ANA;
//======================================================================
typedef VarLoop::Variable Var;
typedef Var2DLoop::Variable2D Var2D;


void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID=1, int targetZ=26, const string playlist="minervame1A", bool doDIS=true);
 
//============================================================================================================================
// GetErrorBands 
//============================================================================================================================

std::map<std::string, std::vector<CVUniverse*> > GetErrorBands(
    PlotUtils::ChainWrapper* chain) {
 
  typedef std::map<std::string, std::vector<CVUniverse*> > SystMap;

  SystMap error_bands;

  // CV
  error_bands[std::string("CV")].push_back(new CVUniverse(chain));

  if(RunCodeWithSystematics) {
  	int n_flux_universes =100 ;
  	SystMap flux_systematics = 
      	PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain,n_flux_universes);
  	error_bands.insert(flux_systematics.begin(), flux_systematics.end());

	//GENIE
	SystMap genie_systematics = 
	     PlotUtils::GetGenieSystematicsMap<CVUniverse>(chain);// change that true to a switch on do_nonrespi_tune
	error_bands.insert(genie_systematics.begin(), genie_systematics.end());
	    
	//Muon Angle systematics...
	SystMap angle_systematics = 
	    PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);//NSFDefaults::beamThetaX_Err,NSFDefaults::beamThetaY_Err);
	error_bands.insert(angle_systematics.begin(), angle_systematics.end());
	  
	//Muon P Systematics...
	//0.01 ---->Why>
	SystMap muonminosP_systematics = PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
	error_bands.insert(muonminosP_systematics.begin(), muonminosP_systematics.end());
	  
	SystMap muonminervaP_systematics = PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
	error_bands.insert(muonminervaP_systematics.begin(), muonminervaP_systematics.end()); 

	//MuonMomentum Resolution Systematics....
	SystMap muonP_resolutions = PlotUtils::GetMuonResolutionSystematicsMap<CVUniverse>(chain);
	error_bands.insert(muonP_resolutions.begin(),muonP_resolutions.end());
	  
	//Minos Efficiency Systematics....
	SystMap minos_efficiency = PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
	error_bands.insert(minos_efficiency.begin(),minos_efficiency.end());

	//2p2h 
	SystMap _2p2h_systematics = PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
	error_bands.insert(_2p2h_systematics.begin(),_2p2h_systematics.end());

	//RPA
	SystMap RPA_systematics = PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
	error_bands.insert(RPA_systematics.begin(),RPA_systematics.end());  
  }

  return error_bands;
}

//============================================================================================================================
// FillVariable 
//============================================================================================================================
    
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
   // std::map< std::string, std::vector<CVUniverse*> > error_bands = utils->GetErrorBands(chain)
   
   std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
   
   std::map<std::string, std::vector<CVUniverse*> >::iterator itr;
   
   std::map<const std::string, const int> error_name;
   for(itr = error_bands.begin(); itr != error_bands.end(); ++itr) error_name.insert(pair<std::string, const int>((itr->second)[0]->ShortName(), (itr->second).size())); 

   std::map<std::string, const int>::iterator itr_m;
   
   std::vector<double> ThetaMuBin, Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin;

   if (doDIS){
     Enubin  = binsDef->GetDISBins("Enu"); 
     Emubin  = binsDef->GetDISBins("Emu"); 
     Ehadbin = binsDef->GetDISBins("Ehad");
     Q2bin = binsDef->GetDISBins("Q2");
     Wbin = binsDef->GetDISBins("W");
     xbin    = binsDef->GetDISBins("x");
     ybin    = binsDef->GetDISBins("y");
     ThetaMuBin = binsDef->GetDISBins("ThetaMu");
     }
   else{
     Enubin  = binsDef->GetEnergyBins("Enu"); 
     Emubin  = binsDef->GetEnergyBins("Emu"); 
     Ehadbin = binsDef->GetEnergyBins("Ehad");
     Q2bin = binsDef->GetEnergyBins("Q2");
     Wbin = binsDef->GetEnergyBins("W");
     xbin    = binsDef->GetEnergyBins("x");
     ybin    = binsDef->GetEnergyBins("y");
   }
  
   //Q2bin = binsDef->GetSidebandBins("Q2");
   //Wbin  = binsDef->GetSidebandBins("W");
   //For 1D varriable

   Var* thetaMu = new Var("GetThetamuDeg", "GetThetamuDeg (Degree)", ThetaMuBin, &CVUniverse::GetThetamuDeg, &CVUniverse::GetThetamuTrueDeg);
   Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
   Var* ehad = new Var("Ehad", "Ehad (GeV)", Ehadbin, &CVUniverse::GetEhadGeV, &CVUniverse::GetEhadTrueGeV);
   Var* Q2 = new Var("Q2", "Q2 (GeV^2)", Q2bin, &CVUniverse::GetQ2RecoGeV, &CVUniverse::GetQ2TrueGeV);
   Var* W = new Var("W", "W (GeV)", Wbin, &CVUniverse::GetWRecoGeV, &CVUniverse::GetWTrueGeV);
   Var* emu = new Var("Emu", "Emu (GeV)", Emubin, &CVUniverse::GetMuonEGeV, &CVUniverse::GetMuonETrueGeV);
   Var* x = new Var("x", "x", xbin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
   Var* y = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);


   
   //std::vector<Var*> variables = {enu,ehad}; 
   variables = { x, y, emu, Q2, enu};//{enu,ehad}; 
   
   //For 2D variable

   //Var2D* W_Q2     = new Var2D(*W, *Q2);
   Var2D* Q2_W     = new Var2D(*Q2, *W);
   Var2D* ThetaMu_emu     = new Var2D(*thetaMu, *emu);
   Var2D* enu_ehad = new Var2D(*enu, *ehad);
   Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
   Var2D* x_Q2 = new Var2D(*x, *Q2);  // y var
   Var2D* y_Q2 = new Var2D(*y, *Q2);  // y var
   Var2D* x_y = new Var2D(*x, *y);  // y var

   variables2d = {emu_ehad, x_y, enu_ehad, Q2_W, y_Q2};//{enu_ehad, Q2_W};
   //variables2d = {emu_ehad};//{enu_ehad, Q2_W};
   
   for (auto v : variables2d) v->InitializeAllHistograms(error_bands);
   for (auto v : variables) v->InitializeAllHistograms(error_bands);

   //Migration starts here!
   for (auto v : variables2d) v->SetupResponse(error_name);
   for (auto v : variables) v->SetupResponse1D(error_name);
    

   int mc0=0;
   int mc1=0;
   int mc2=0;
   int mc3=0;
   int mc4=0; 
   int mc5=0; 
   int mc6=0; 
   int mc7=0; 
   int mc8=0; 
   int mc9=0; 
   int mc10=0; 
   int mc11=0; 
   int mc12=0; 
   int mc13=0; 
   int mc14=0; 
   int mc15=0; 
   int mc16=0; 
   int mc17=0; 
   
   CVUniverse *dataverse = new CVUniverse(chain,0);
    
   //=========================================
   // Entry Loop
   //=========================================

   std::cout<<"# of entries = "<<chain->GetEntries()<<std::endl;

   for(int i=0; i<chain->GetEntries(); ++i){
     if(i%500000==0) std::cout << (i/1000) << "k " << std::endl;
     //=========================================
     // For every systematic, loop over the universes, and fill the
     // appropriate histogram in the MnvH1D
     //=========================================

       for (auto band : error_bands){
	       int unv_count = 0;
               std::vector<CVUniverse*> error_band_universes = band.second;
	       for (auto universe : error_band_universes){
		 // Tell the Event which entry in the TChain it's looking at
		 //=========================================
		 // CUTS in each universe
		 //========================================
	 
		 universe->SetEntry(i);
	   if(!cutter->PassHelicityCut(universe,helicity)) continue;
           mc0++;
	   if( ! universe->GetInt("NukeCC_in_fiducial_area") ) continue;
           mc1++;
	   if( ! cutter->PassZDistCut( universe) ) continue;
           mc2++;
	   if( ! cutter->PassDistToDivisionCut(universe) ) continue;
           mc3++;
	   if( 1 < universe->GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj") ) continue;
           mc4++;
	   if( 120. < universe->GetEnu() * mev_to_gev ) continue;
           mc5++;
	   if( ! cutter->PassMuCurveCut(universe, helicity) ) continue;
           mc6++;
	   if( ! cutter->PassMuCoilCut(universe) ) continue;
           mc7++;
	   if(!cutter->PassReco(universe,helicity)) continue;
           mc8++;
	   if(!cutter->IsInMaterial(universe,targetID,targetZ, false)) continue;
           mc9++;
	   if(targetID<10 && universe->GetInt("NukeCC_targetID") != targetID) continue;
           mc10++;
	   if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
	   mc11++;
           if(!cutter->PassDISCut(universe)) continue;
	   
           if (!cutter->PassTrueFiducial(universe)) continue;
           mc12++;
	   if(!cutter->IsInTrueMaterial(universe,targetID, targetZ,false)) continue;
           mc13++;
           if(!cutter->PassTrueDISCut(universe)) continue;
		 
          for (auto v : variables2d){
	     if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(universe)) continue;
	     if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(universe)) continue;
	     if( v->GetNameX()=="Emu")mc14++;
             if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassTrueMuEnergyCut(universe)) continue;
	     if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassTrueThetaCut(universe)) continue;
	     //double wgt = universe->GetWeight();
             if( v->GetNameX()=="Emu")mc15++;
             v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeighty()); 
             //v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), wgt); 
             v->m_selected_Migration.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetTrueValueX(*universe), universe->GetWeighty()); 
             //v->m_selected_Migration.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), wgt); 
     //Migration stuff
             v->FillResponse(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),universe->ShortName(),universe->GetWeighty(),unv_count); 
		 }

	   for (auto v : variables){
	     if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(universe)) continue;
	     if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(universe))continue;
	     if( v->GetName()=="Emu")mc16++;
             if( v->GetName()!="Emu")   if(!cutter->PassTrueMuEnergyCut(universe)) continue;
	     if( v->GetName()!="ThetaMu") if(!cutter->PassTrueThetaCut(universe))continue;
	     if( v->GetName()=="Emu")mc17++;
	     //Added for 1D analysis 
             v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeighty());
	     v->m_selected_mc_true.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetWeighty());
	    //////////////Migration//////////////////
              v->FillResponse1D(v->GetRecoValue(*universe),v->GetTrueValue(*universe),universe->ShortName(),universe->GetWeighty(),unv_count); 
	   	 }
		 unv_count++;
               } // End band's universe loop
     }
  }//End entries loop
 
 for (auto v : variables2d) v->getResponseObjects(error_bands);
 for (auto v : variables) v->getResponseObjects1D(error_bands);
 
 for(auto band : error_bands){
    std::vector<CVUniverse*> band_universes = band.second;
    for(unsigned int i_universe = 0; i_universe < band_universes.size(); ++i_universe)
      delete band_universes[i_universe];
 } 
    
 delete dataverse;

   std::cout<<"**********************************"<<std::endl;
      std::cout<<" Summary "<<std::endl;
   std::cout<<"**********************************"<<std::endl;
     std::cout<<"Printing the Numerator Summary "<<std::endl;
     std::cout<<" Reco cut0 = "<<mc0<<std::endl;
     std::cout<<" Reco cut1 = "<<mc1<<std::endl;
     std::cout<<" Reco cut2 = "<<mc2<<std::endl;
     std::cout<<" Reco cut3 = "<<mc3<<std::endl;
     std::cout<<" Reco cut4 = "<<mc4<<std::endl;
     std::cout<<" Reco cut5 = "<<mc5<<std::endl;
     std::cout<<" Reco cut6 = "<<mc6<<std::endl;
     std::cout<<" Reco cut7 = "<<mc7<<std::endl;
     std::cout<<" Reco cut8 = "<<mc8<<std::endl;
     std::cout<<" Reco cut9 = "<<mc9<<std::endl;
     std::cout<<" Reco cut10 = "<<mc10<<std::endl;
     std::cout<<" Reco cut11= "<<mc11<<std::endl;
     std::cout<<" Reco cut12 = "<<mc12<<std::endl;
     std::cout<<" Reco cut13= "<<mc13<<std::endl;
     std::cout<<" Reco cut14= "<<mc14<<std::endl;
     std::cout<<" Reco cut15= "<<mc15<<std::endl;
     std::cout<<" Reco cut16= "<<mc16<<std::endl;
     std::cout<<" Reco cut17= "<<mc17<<std::endl;
}

//============================================================================================================================
// Main
//============================================================================================================================

int main(int argc, char *argv[]){
	 ROOT::Cintex::Cintex::Enable();
	 TH1::AddDirectory(false);

	 if(argc==1){
	     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
	     std::cout<<"MACROS HELP:\n\n"<<
	       "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n"<<
	       "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
	       "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
	       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" << std::endl;
	     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
	     return 0;
	 }

	 TString dir(argv[1]);
	 int targetID = atoi(argv[2]);
	 int targetZ  = atoi(argv[3]);

		
	 // TString dir(argv[1]);
	 // int targetID = 1;
	 // int targetZ = 26;
	 // const string playlist= argv[4];
	   
	 bool doDIS = true;
	 bool isSYS = true;
	 const std::string mc_file_list("../../include/playlists/NukeCC_minervame1A_MC_Inextinguishable_merged.txt");
	 const std::string data_file_list("../../include/playlists/NukeCC_minervame1A_DATA_Inextinguishable_merged.txt");
	 const std::string plist_string("minervame1A");
	 const std::string reco_tree_name("NukeCC");
	 const bool wants_truth = false;
	 const bool is_grid = false;

	 // PlotUtils::MacroUtil util(reco_tree_name, file_list, plist_string, wants_truth,
	 //                         is_grid);
	 PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth, is_grid);

	 util.PrintMacroConfiguration("main");

	 //=========================================
	 // Systematics
	 //=========================================
	 //std::map<std::string, std::vector<CVUniverse*> > error_bands =
	 //  GetErrorBands(util.m_mc);

	 PlotUtils::DefaultCVUniverse::SetNFluxUniverses(100);
	 PlotUtils::DefaultCVUniverse::SetNuEConstraint(true);
	 PlotUtils::DefaultCVUniverse::SetAnalysisNuPDG(14);
	 PlotUtils::DefaultCVUniverse::SetNonResPiReweight(true);
	    
	 NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
	 NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
	 NukeCC_Binning  *binsDef = new NukeCC_Binning();
	  
	 PlotUtils::ChainWrapper* chainData = util.m_data;
	 PlotUtils::ChainWrapper* chainMC = util.m_mc;
	 HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);
	 double DataPot=  util.m_data_pot; 
	 double MCPot=  util.m_mc_pot; 
	 //double total_pot_data,total_pot_mc;
	 //utils->getPOT(total_pot_data,total_pot_mc);  
	 double  MCscale=DataPot/MCPot;
	 //double  MCscale=1.0;
	 
	 std::cout<<" MCScale= "<<MCscale<<std::endl; 
	 std::vector<Var*> variablesMC,variablesData; 
	 std::vector<Var2D*> variables2DMC,variables2DData; 

	 TString histFileName = utils->GetHistFileName( "Migrationy", FileType::kAny, targetID, targetZ, helicity ); 
	   
	 TFile fout(dir.Append(histFileName),"RECREATE");	
	   
	 // For 1D variables 
	 FillVariable(chainMC, helicity, utils, cutter, binsDef, variablesMC, variables2DMC, true, targetID, targetZ, plist_string, doDIS);
	     
         	 
	 for (auto v : variablesMC) v-> mresp1D.SyncCVHistos();
	 for (auto v : variables2DMC) v-> mresp.SyncCVHistos();
	 
	 
	 for (auto v : variablesMC) {
	   v->WriteAllHistogramsToFile(fout, true);
	 }


	 // Plotting If you want for 1D
	 /* 
	 for(int i=0; i< variablesMC.size();i++){
		 PlotCVAndError(variablesData[i]->m_selected_data_reco.hist,variablesMC[i]->m_selected_mc_reco.hist,variablesMC[i]->GetName(),MCscale);
		       
		 PlotErrorSummary(variablesMC[i]->m_selected_mc_reco.hist, variablesMC[i]->GetName());
		 PlotStacked(variablesData[i]->m_selected_data_reco_sb.hist,variablesMC[i]->m_selected_mc_sb.GetHistArray(),MCscale, variablesMC[i]->m_selected_mc_sb.GetName(), variablesMC[i]->m_selected_mc_sb.GetName());
	 }//End 1D plotting 
	 */

	 //For 2D variable

	 for (auto v : variables2DMC) {
	   v->WriteAllHistogramsToFile(fout,true);
	 }

	 
	 //Plotting in 2D
	   
	 for(int i=0; i< variables2DMC.size();i++){
	   Plot2D(variables2DMC[i]->mresp.hist, variables2DMC[i]->GetName(), variables2DMC[i]->GetNameX(), variables2DMC[i]->GetNameY()); //Plotting line that I somehow cannot delete without producing memory errors, but no one else can reproduce. --ANF 2020.4.6
	 //Plot2D(variables2DData[i]->m_selected_data_reco.hist, variables2DData[i]->GetName(), variables2DData[i]->GetNameX(),variables2DData[i]->GetNameY());
	     
	 }//End 2D plotting

}//End Main

