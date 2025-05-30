//==============================================================================
// Loop entries, make cuts, fill histograms.
// * Uses the New Systematics Framework and "Universe" objects.
// * loop universes, make cuts and fill histograms with the correct lateral
// shifts and weights for each universe.
// * TChain --> PlotUtils::ChainWrapper.
// * MnvHXD --> PlotUtils::HistWrapper.
// * Genie, flux, non-resonant pion, and some detector systematics calculated.
//==============================================================================

//#include "../include/Variable.h"
#include "../../../NUKECCSRC/ana_common/include/CommonIncludes.h"
#include "../../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "../../include/Variable.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"   
#include "../../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include <iostream>
#include <stdlib.h>
#include "../../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../../NUKECCSRC/ana_common/include/NukeCC_Cuts.h"
#include "TParameter.h"

#include "../../include/systematics/Systematics.h"

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
// Main
//============================================================================================================================

int main(int argc, char *argv[]){
	 //ROOT::Cintex::Cintex::Enable();
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
   const string playlist= argv[4];


		bool doDIS = false; 
	 // TString dir(argv[1]);
	 // int targetID = 1;
	 // int targetZ = 26;
	 // const string playlist= argv[4];

    const std::string plist_string(playlist);
    const std::string mc_file_list(Form("../../include/playlists/MasterAnaDev_MC_%s_DSCAL_MuonKludged.txt", plist_string.c_str()));
    const std::string data_file_list(Form("../../include/playlists/MasterAnaDev_Data_%s_DSCAL_MuonKludged.txt",plist_string.c_str()));
    const std::string reco_tree_name("MasterAnaDev");
  
    const bool wants_truth = false;
    //const bool is_grid = false;
    // is grid removed after update of MAT 07/12/2021

	 PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth);

	 util.PrintMacroConfiguration("main");

	 //=========================================
	 // Systematics
	 //=========================================
	 //std::map<std::string, std::vector<CVUniverse*> > error_bands =
	 //  GetErrorBands(util.m_mc);

  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
  PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);
  // Defined for MnvHadronReweighter (GEANT Hadron sytematics)
  //Tracker or nuke (what clusters are accepted for reconstruction)
  PlotUtils::MinervaUniverse::SetReadoutVolume("Tracker");
  //Neutron CV reweight is on by default (recommended you keep this on)
  //PlotUtils::MinervaUniverse::SetMHRWeightNeutronCVReweight(true);
  //Elastics are on by default (recommended you keep this on)
  //PlotUtils::MinervaUniverse::SetMHRWeightElastics(true);
  PlotUtils::MinervaUniverse::RPAMaterials(true);
	    
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

         //TString histFileName = utils->GetHistFileName( "Migration", FileType::kAny, targetID, targetZ );
	TString histFileName;
  if(RunCodeWithSystematics){
    histFileName += Form("/Migration_%s_t%d_z%02d_sys.root", plist_string.c_str(), targetID, targetZ);
  }

  else{
    histFileName += Form("/Migration_%s_t%d_z%02d_nosys.root", plist_string.c_str(), targetID, targetZ);
  } 
  	   
	TFile fout(dir.Append(histFileName),"RECREATE");	
	   
	 // For 1D variables 
	 FillVariable(chainMC, helicity, utils, cutter, binsDef, variablesMC, variables2DMC, true, targetID, targetZ, plist_string, doDIS);
	     
         	 
	 for (auto v : variablesMC) v-> mresp1D.SyncCVHistos();
   for (auto v : variablesMC) v->m_selected_mc_reco.SyncCVHistos();   
   for (auto v : variablesMC) v->m_selected_Migration.SyncCVHistos(); 

	 for (auto v : variables2DMC) v-> mresp.SyncCVHistos();
	 
	 
	 for (auto v : variablesMC) {
	  v->WriteAllHistogramsToFileMig(fout, true);  
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
	   v->WriteAllHistogramsToFileMig(fout,true);
	 }

  fout.cd();
  auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  dataPOTOut->Write();
  mcPOTOut->Write(); 
	 
	 //Plotting in 2D
	   
	 //for(int i=0; i< variables2DMC.size();i++){
	   //Plot2D(variables2DMC[i]->mresp.hist, variables2DMC[i]->GetName(), variables2DMC[i]->GetNameX(), variables2DMC[i]->GetNameY()); //Plotting line that I somehow cannot delete without producing memory errors, but no one else can reproduce. --ANF 2020.4.6
	 //Plot2D(variables2DData[i]->m_selected_data_reco.hist, variables2DData[i]->GetName(), variables2DData[i]->GetNameX(),variables2DData[i]->GetNameY());
	     
	 //}//End 2D plotting
  std::cout << "DONE" << std::endl;

}//End Main

//============================================================================================================================
// FillVariable 
//============================================================================================================================
    
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
   // std::map< std::string, std::vector<CVUniverse*> > error_bands = utils->GetErrorBands(chain)
   
   std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
   
   std::map<std::string, std::vector<CVUniverse*> >::iterator itr;
   
   std::map<const std::string, int> error_name;
   for(itr = error_bands.begin(); itr != error_bands.end(); ++itr) error_name.insert(pair<std::string, const int>((itr->second)[0]->ShortName(), (itr->second).size())); 

   std::map<std::string, const int>::iterator itr_m;



//===========================DEBUG TO KNOW ABOUT THE UNIVERSES=========================================================
if (error_bands.empty()) {
    std::cerr << "[ERROR] error_bands is EMPTY!" << std::endl;
} else {
    std::cout << "[INFO] error_bands contains " << error_bands.size() << " systematic error bands." << std::endl;
    std::cout << "[INFO] error_name contains " << error_name.size() << " systematic error band names." << std::endl;
}

std::cout << "\n[INFO] Contents of error_name:\n";
for (const auto& [short_name, n_universes] : error_name) {
    std::cout << "  Error band: " << short_name << " - Number of universes: " << n_universes << std::endl;
}


for (const auto& [name, universes] : error_bands) {
    std::cout << "Systematic: " << name << " - Number of universes: " << universes.size() << std::endl;

    if (universes.empty()) {
        std::cerr << "[WARNING] Systematic " << name << " has no universes!" << std::endl;
    } else if (!universes[0]) {
        std::cerr << "[WARNING] First universe in " << name << " is a nullptr!" << std::endl;
    } else {
        std::cout << "  First universe ShortName: " << universes[0]->ShortName() << std::endl;
    }
}
//======================================================================================================================

   std::vector<double> ThetaMuBin, Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin, xbinBrian,Eavbin,q3bin;
   std::vector<double> x09bin, xfinebin;
   std::vector<double> pTbin, pZbin, pZbin1D, pTbin1D;

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
     Eavbin=binsDef->GetEnergyBins("Eavail");
     q3bin =binsDef->GetEnergyBins("q3");
     x09bin = binsDef->GetEnergyBins("x09");
     xfinebin = binsDef->GetEnergyBins("xfine");
     xbinBrian    = binsDef->GetEnergyBins("xBrian");
     ybin    = binsDef->GetEnergyBins("y");
     pTbin = binsDef->GetEnergyBins("muonPt"); 
     pZbin = binsDef->GetEnergyBins("muonPz"); 
     pZbin1D = binsDef->GetEnergyBins("muonPz1D"); 
     pTbin1D = binsDef->GetEnergyBins("muonPt1D");
     ThetaMuBin = binsDef->GetEnergyBins("ThetaMu");
   }
  
   //Q2bin = binsDef->GetSidebandBins("Q2");
   //Wbin  = binsDef->GetSidebandBins("W");
   //For 1D varriable

   //Var* thetaMu = new Var("ThetamuDeg", "ThetamuDeg", ThetaMuBin, &CVUniverse::GetThetamuDeg, &CVUniverse::GetThetamuTrueDeg);
   //Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
   //Var* ehad = new Var("Ehad", "Ehad (GeV)", Ehadbin, &CVUniverse::GetEhadGeV, &CVUniverse::GetEhadTrueGeV);
   //Var* Q2 = new Var("Q2", "Q2 (GeV^2)", Q2bin, &CVUniverse::GetQ2RecoGeV, &CVUniverse::GetQ2TrueGeV);
   //Var* W = new Var("W", "W (GeV)", Wbin, &CVUniverse::GetWRecoGeV, &CVUniverse::GetWTrueGeV);
   //Var* emu = new Var("Emu", "Emu (GeV)", Emubin, &CVUniverse::GetMuonEGeV, &CVUniverse::GetMuonETrueGeV);
   //Var* x = new Var("x", "x", xbin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
   //Var* x09 = new Var("x09", "x09", x09bin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
   //Var* xfine = new Var("xfine", "xfine", xfinebin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
   //Var* xBrian = new Var("xBrian", "xBrian", xbinBrian, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
   //Var* y = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);
   //Var* pTmu = new Var("pTmu", "pTmu", pTbin, &CVUniverse::GetMuonPt, &CVUniverse::GetlepPtTrue);
   //Var* pZmu = new Var("pZmu", "pZmu", pZbin, &CVUniverse::GetMuonPz, &CVUniverse::GetlepPzTrue);
   //Var* pTmu1D = new Var("pTmu1D", "pTmu1D", pTbin1D, &CVUniverse::GetMuonPt, &CVUniverse::GetlepPtTrue);
   //Var* pZmu1D = new Var("pZmu1D", "pZmu1D", pZbin1D, &CVUniverse::GetMuonPz, &CVUniverse::GetlepPzTrue);
   
   Var* Eavail = new Var("Eavailable", "Eavail (GeV)", Eavbin, &CVUniverse::GetTrackerECALAvEnergy, &CVUniverse::GetEavailTrue);
   Var* q3 = new Var("q3", "q3 (GeV)", q3bin, &CVUniverse::Getq3Reco, &CVUniverse::Getq3Truth);

   //std::vector<Var*> variables = {enu,ehad}; 
   variables = {Eavail, q3};//{enu,ehad}; 
   
   //For 2D variable

    //Var2D* W_Q2     = new Var2D(*W, *Q2);
    //Var2D* enu_ehad = new Var2D(*enu, *ehad);
    //Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
    //Var2D* x_Q2 = new Var2D(*x, *Q2);  // y var
    //Var2D* x_y = new Var2D(*x, *y);  // y var
   //Var2D* pZmu_pTmu = new Var2D(*pZmu, *pTmu);
    Var2D* Eav_q3=new Var2D(*Eavail, *q3);
    
    variables2d = {Eav_q3};//{pZmu_pTmu };
   
   for (auto v : variables2d) v->InitializeAllHistograms(error_bands);
   for (auto v : variables) v->InitializeAllHistograms(error_bands);

   //Migration starts here!
   for (auto v : variables2d) v->SetupResponse(error_name);
   for (auto v : variables) v->SetupResponse1D(error_name);
    
  int reco0=0;
  int reco1=0;
  int reco2=0;
  int reco3=0;
  int reco4=0; 
  int reco5=0; 
  int reco6=0;  
  int reco7 = 0;
  int reco8=0;
  

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
     reco0++;
     // reco cuts
     if( universe->GetInt("muon_corrected_p") == -999 ) continue; // additional cut to get rid of an issue
     reco1++;
     if(!cutter->PassReco(universe,helicity)) continue;
     reco2++;
     int rock_muons_removed = universe->GetInt("rock_muons_removed");
     if(rock_muons_removed == 1) continue;
     //if(!cutter->IsInMaterial(universe,targetID,targetZ, false)) continue;
     reco3++;
     if(!cutter->isFiducialRegion(universe,5990.0 , 8340.0, 850.0)) continue; // cv, zmin,zmax,apothem
     //if(targetID<10 && universe->GetInt("MasterAnaDev_ANN_targetID") != targetID) continue;
     reco4++;
     if (!cutter->PassTrueCC(universe, helicity)) continue; //true CC, true nu
     //if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
     reco5++;
     if(!cutter->PassMuEnergyCut(universe)) continue;
      // NO xy separation,  NO APOTHEM CUT
      reco6++;

	   //if(!cutter->IsInTrueMaterial(universe,targetID, targetZ,false)) continue; // true target + material
     reco7++;

		 for (auto v : variables2d){
	    //if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(universe)) continue;
	    //if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(universe)) continue;
      
      //if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassTrueMuEnergyCut(universe)) continue;
	    //if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassTrueThetaCut(universe)) continue;
      // NO TRUE angle cut, efficiency corrected
	  
      
      v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
      v->m_selected_Migration.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight()); 
      
      //Migration stuff
      v->FillResponse(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),universe->ShortName(),universe->GetWeight(),unv_count); 
		 }

	   for (auto v : variables){
	     //if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(universe)) continue;
	     //if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(universe))continue;
	     //if( v->GetName()=="Enu") reco7++;
       
       //if( v->GetName()!="Emu")   if(!cutter->PassTrueMuEnergyCut(universe)) continue;
	     //if( v->GetName()!="ThetaMu") if(!cutter->PassTrueThetaCut(universe))continue;
       // NO TRUE angle cut, efficiency corrected
	     //if( v->GetName()=="Enu") reco8++;

	     v->m_selected_mc_reco.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetWeight());
       v->m_selected_Migration.univHist(universe)->Fill(v->GetRecoValue(*universe), v->GetTrueValue(*universe), universe->GetWeight()); 
       //1D response
       v->FillResponse1D(v->GetRecoValue(*universe),v->GetTrueValue(*universe),universe->ShortName(),universe->GetWeight(),unv_count); 
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
     std::cout<<"Migration Matrix "<<std::endl;
    std::cout << "No cuts = " << reco0 << std::endl;
    std::cout << "Reco Cut = " << reco1 << std::endl;
    std::cout << "Material Cut = " << reco2 << std::endl;
    std::cout << "TargetID Cuts = " << reco3 << std::endl;
    std::cout << "Plane prob. cut = " << reco4 << std::endl;
    std::cout << "Truth cut (fiducial, CC, antinu) = " << reco5 << std::endl;
    std::cout << "True Material cut  = "<< reco6 << std::endl;
    std::cout << "Muon kinematics cut  = " << reco7 << std::endl;
    std::cout << "True muon kinematics cut  = " << reco8 << std::endl;
}


