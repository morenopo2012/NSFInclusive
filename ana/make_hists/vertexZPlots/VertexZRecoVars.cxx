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
#include "../../../NUKECCSRC/ana_common/include/CommonIncludes.h"
#include "../../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "Variable.h"
#include "../../include/playlists/playlists.h"
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
#include "../../../NUKECCSRC/ana_common/include/NukeCC_Cuts.h"
#include "TParameter.h"


#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
//#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
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

int main(int argc, char *argv[]){
    //ROOT::Cintex::Cintex::Enable();
    TH1::AddDirectory(false);

 if(argc==1){
   std::cout<<"---------------------------------------------------------------------------------"<<std::endl;
   std::cout<<"MACROS HELP:\n\n"<<
     "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n"<<
     "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
     "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
     "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" << std::endl;
     std::cout<<"-------------------------------------------------------------------------------"<<std::endl;
     return 0;
   }

  TString dir(argv[1]);
  int targetID = 99;
  int targetZ =  99;

  const string plist_string= argv[4];
  bool doDIS=false;
  const std::string mc_file_list( get_mc_files(plist_string, targetID) );
  const std::string data_file_list( get_data_files(plist_string) );
 
  //Keep this here it is for older way of caling anatuples
  //const std::string data_file_list("../../include/playlists/MasterAnaDev_Data_minervame6B_MuonKludged.txt");
  //const std::string mc_file_list("../../include/playlists/MasterAnaDev_MC_minervame6B_MuonKludged.txt");

  //const std::string plist_string("minervame6B");
  const std::string reco_tree_name("MasterAnaDev");
  const bool wants_truth = false;
  const bool is_grid = false;

   PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth);
   //PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, playlist, wants_truth);

   util.PrintMacroConfiguration("main");

  //=========================================
  // Systematics
  //=========================================
  //std::map<std::string, std::vector<CVUniverse*> > error_bands =
  //GetErrorBands(util.m_mc);

   PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
   PlotUtils::MinervaUniverse::SetNuEConstraint(true);
   PlotUtils::MinervaUniverse::SetAnalysisNuPDG(-14);
   PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
   PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
   PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);
   PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);

     
   NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
   NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
   NukeCC_Binning  *binsDef = new NukeCC_Binning();
  
   PlotUtils::ChainWrapper* chainData = util.m_data;
   PlotUtils::ChainWrapper* chainMC = util.m_mc;
   HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);
   double DataPot=  util.m_data_pot; 
   double MCPot=  util.m_mc_pot; 
   double  MCscale=DataPot/MCPot;

 
   std::cout<<" MCScale= "<<MCscale<<std::endl; 
   std::vector<Var*> variablesMC,variablesData; 
   std::vector<Var2D*> variables2DMC,variables2DData; 
  

   TString histFileName = utils->GetHistFileNamePlaylist( "EventSelection", FileType::kAny, targetID, targetZ, helicity, plist_string );
  
   //Works good for the grid submission
   //TString histFileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ );
   
   TFile fout(dir.Append(histFileName),"RECREATE");	
   
   // For 1D variables 
   FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true,targetID, targetZ, plist_string,doDIS);
       
   for (auto v : variablesMC) {
    v->m_selected_mc_reco_water.SyncCVHistos();
    v->m_selected_mc_reco_carbon.SyncCVHistos();
    v->m_selected_mc_reco_lead.SyncCVHistos();  
    v->m_selected_mc_reco_iron.SyncCVHistos();
    v->m_selected_mc_reco_scintillator.SyncCVHistos();

  }
   
   FillVariable(chainData, helicity, utils, cutter,binsDef,variablesData,variables2DData,false,targetID, targetZ, plist_string,doDIS);
   
   for (auto v : variablesData) {
    v->m_selected_data_reco.SyncCVHistos();
    v->m_selected_data_reco_water.SyncCVHistos();
    v->m_selected_data_reco_carbon.SyncCVHistos();
    v->m_selected_data_reco_lead.SyncCVHistos();  
    v->m_selected_data_reco_iron.SyncCVHistos();
    v->m_selected_data_reco_scintillator.SyncCVHistos();
    }

   //for (auto v : variables2DData) v->m_selected_data_reco.SyncCVHistos();
 
   for (auto v : variablesMC) {
     v->WriteAllHistogramsToFile(fout, true);
   }


   for (auto v : variablesData) {
     v->WriteAllHistogramsToFile(fout, false);
   }


  //Writing POT to the HistFile
  fout.cd();
  auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  dataPOTOut->Write();
  mcPOTOut->Write(); 
  
}//End Main

   
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
  std::map<std::string, std::vector<CVUniverse*> > error_bands =GetErrorBands(chain);
  std::vector<double> vtxzbin;


  if (doDIS){
     //vtxzbin = binsDef->GetEnergyBins("vtxz_all");
     }
  else{
     vtxzbin = binsDef->GetEnergyBins("vtxz_all");
   }

   
   //Var* vtxz = new Var("vtxz", "Vertex Z (cm)", vtxzbin, &CVUniverse::GetVertexZMy, &CVUniverse::GetVertexZTrueMy);
   Var* vtxz = new Var("vtxz", "Vertex Z (cm)", vtxzbin, &CVUniverse::GetVertexZNew, &CVUniverse::GetVertexZTrueMy);



   variables = {vtxz};//{enu,ehad}; 
   
   
   for (auto v : variables) v->InitializeAllHistograms(error_bands);
   //for (auto v : variables2d) v->InitializeAllHistograms(error_bands);


   double no_cuts = 0;
   int reco0=0;
   int reco1=0;
   int reco2=0;
   int reco3=0;
   int reco4=0; 
  int scintillator=0;
  int water=0;
  int iron=0;
  int lead=0;
  int carbon=0; 

   
   CVUniverse *dataverse = new CVUniverse(chain,0);
    
 //  bool anyTrackerMod = false; 
   //=========================================
   // Targets combining Loop
   //=========================================
  ////for(int t = 0; t < targetIDs.size(); t++){	  		
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
       if(isMC){     
       for (auto band : error_bands){
       std::vector<CVUniverse*> error_band_universes = band.second;
       for (auto universe : error_band_universes){
	 // Tell the Event which entry in the TChain it's looking at
	 //=========================================
	 // CUTS in each universe
	 //========================================
	 universe->SetEntry(i);
          reco0++;
            if( ! cutter->PassTrueCC(universe,helicity)) continue;
            //if( ! cutter->InHexagonTrue(universe, 850.) ) continue;
            if( ! cutter->InHexagon(universe, 850.) ) continue;
            //if(!cutter->PassReco(universe,helicity)) continue;   
            reco1++;

            //if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
            //if( universe->GetVecElem("ANN_plane_probs",0) < 0.2 ) continue;	   
            reco2++;

            if(!cutter->PassMuEnergyCut(universe)) continue;
            reco3++;

            if(!cutter->PassThetaCut(universe))continue;
            reco4++;

            // WATER
            if(cutter->Water(universe)){
            //if(cutter->WaterTrue(universe)){
            //if(universe->GetInt("MasterAnaDev_targetID") ==2){
              water++;
              for (auto v : variables){
              v->m_selected_mc_reco_water.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              }
            }
            
            // CARBON
            //else if(cutter->IsInTrueMaterial(universe,3,6, /*anyTrakerMod*/false)){
            else if(cutter->IsInMaterial(universe,3,6, /*anyTrakerMod*/false)){
              carbon++;
              for (auto v : variables){
                v->m_selected_mc_reco_carbon.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              }
            }

            //IRON T1235
            //else if( (cutter->IsInTrueMaterial(universe,1,26, /*anyTrakerMod*/false)) || (cutter->IsInTrueMaterial(universe,2,26, /*anyTrakerMod*/false))
            //    || (cutter->IsInTrueMaterial(universe,3,26, /*anyTrakerMod*/false)) ||  (cutter->IsInTrueMaterial(universe,5,26, /*anyTrakerMod*/false))){
            else if( (cutter->IsInMaterial(universe,1,26, /*anyTrakerMod*/false)) || (cutter->IsInMaterial(universe,2,26, /*anyTrakerMod*/false))
                || (cutter->IsInMaterial(universe,3,26, /*anyTrakerMod*/false)) ||  (cutter->IsInMaterial(universe,5,26, /*anyTrakerMod*/false))){
            //else if( (cutter->IsInTrueMaterial(universe,3,26, /*anyTrakerMod*/false))){
              iron++;
              for (auto v : variables){
                v->m_selected_mc_reco_iron.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              }
            }
         

            // LEAD T12345
            //else if( (cutter->IsInTrueMaterial(universe,1,82, /*anyTrakerMod*/false)) || (cutter->IsInTrueMaterial(universe,2,82, /*anyTrakerMod*/false)) 
             //   || (cutter->IsInTrueMaterial(universe,3,82, /*anyTrakerMod*/false)) || (cutter->IsInTrueMaterial(universe,4,82, /*anyTrakerMod*/false)) 
             //   ||  (cutter->IsInTrueMaterial(universe,5,82, /*anyTrakerMod*/false)) ){
            else if( (cutter->IsInMaterial(universe,1,82, /*anyTrakerMod*/false)) || (cutter->IsInMaterial(universe,2,82, /*anyTrakerMod*/false)) 
                || (cutter->IsInMaterial(universe,3,82, /*anyTrakerMod*/false)) || (cutter->IsInMaterial(universe,4,82, /*anyTrakerMod*/false)) 
                ||  (cutter->IsInMaterial(universe,5,82, /*anyTrakerMod*/false)) ){
            //else if((cutter->IsInTrueMaterial(universe,3,82, /*anyTrakerMod*/false))  ){
              lead++;
              for (auto v : variables){
                v->m_selected_mc_reco_lead.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());   
              }
            }

            else if(cutter->TrackerOnly(universe)){
              for (auto v : variables){
                scintillator++;
                v->m_selected_mc_reco_scintillator.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              }
            }

          } // End band's universe loops
        }// End Band loop
      }



        else{

         dataverse->SetEntry(i);
         reco0++;
      if(!cutter->PassReco(dataverse,helicity)) continue;
      reco1++;
      
      //if( dataverse->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
      //if( dataverse->GetVecElem("ANN_plane_probs",0) < 0.2 ) continue;	    
      reco2++;

      if(!cutter->PassMuEnergyCut(dataverse)) continue;   
      reco3++;

      if(!cutter->PassThetaCut(dataverse))continue;
      reco4++;
    

      // WATER
      if(cutter->Water(dataverse)){
      //if(universe->GetInt("MasterAnaDev_targetID") ==2){
        water++;
        for (auto v : variables){
        v->m_selected_data_reco_water.hist->Fill(v->GetRecoValue(*dataverse, 0));	
        }
      }
      
      // CARBON
      else if( (cutter->IsInMaterial(dataverse,3,6, /*anyTrakerMod*/false) && dataverse->GetInt("MasterAnaDev_ANN_targetID") == 3)){
        carbon++;
        for (auto v : variables){
          v->m_selected_data_reco_carbon.hist->Fill(v->GetRecoValue(*dataverse, 0));	

        }
      }
    
      //IRON Ts235
      else if( (cutter->IsInMaterial(dataverse,2,26, /*anyTrakerMod*/false)  && dataverse->GetInt("MasterAnaDev_ANN_targetID") == 2)
          || (cutter->IsInMaterial(dataverse,3,26, /*anyTrakerMod*/false) && dataverse->GetInt("MasterAnaDev_ANN_targetID") == 3) ||  (cutter->IsInMaterial(dataverse,5,26, /*anyTrakerMod*/false) && dataverse->GetInt("MasterAnaDev_ANN_targetID") == 5)){
      //else if( (cutter->IsInTrueMaterial(universe,3,26, /*anyTrakerMod*/false))){
        iron++;
        for (auto v : variables){
          v->m_selected_data_reco_iron.hist->Fill(v->GetRecoValue(*dataverse, 0));	

        }
      }
      
      // LEAD T2345
      else if( (cutter->IsInMaterial(dataverse,2,82, /*anyTrakerMod*/false) && dataverse->GetInt("MasterAnaDev_ANN_targetID") == 2) 
          || (cutter->IsInMaterial(dataverse,3,82, /*anyTrakerMod*/false)  && dataverse->GetInt("MasterAnaDev_ANN_targetID") == 3) || (cutter->IsInMaterial(dataverse,4,82, /*anyTrakerMod*/false) && dataverse->GetInt("MasterAnaDev_ANN_targetID") == 4) 
          ||  (cutter->IsInMaterial(dataverse,5,82, /*anyTrakerMod*/false)  && dataverse->GetInt("MasterAnaDev_ANN_targetID") == 5) ){
      //else if((cutter->IsInTrueMaterial(universe,3,82, /*anyTrakerMod*/false))  ){
        lead++;
        for (auto v : variables){
          v->m_selected_data_reco_lead.hist->Fill(v->GetRecoValue(*dataverse, 0));	

        }
      }

      else if(cutter->TrackerOnly(dataverse)){
        for (auto v : variables){
          scintillator++;
          v->m_selected_data_reco_scintillator.hist->Fill(v->GetRecoValue(*dataverse, 0));	

        }
      }
          
    }

  }//End entries loop


for(auto band : error_bands){
    std::vector<CVUniverse*> band_universes = band.second;
    for(unsigned int i_universe = 0; i_universe < band_universes.size(); ++i_universe)
    delete band_universes[i_universe];
  } 
    
    delete dataverse;


  // std::cout<<"**********************************"<<std::endl;
   std::cout<<"Printing the ";
   isMC? std::cout<<"MC": std::cout<<"Data";
   std::cout<<" Summary "<<std::endl;
   std::cout<<" No cuts = "<<reco0<<std::endl;
   std::cout<<" Reco Cut = "<<reco1<<std::endl;
   std::cout<<" Material Cut = "<<reco2<<std::endl;
   std::cout<<" Plane/targetID Cuts = "<<reco3<<std::endl;
   std::cout<<" Muon Kinematics Cuts = "<<reco4<<std::endl;
   std::cout<<"**********************************"<<std::endl;
   std::cout<<" Summary "<<std::endl;
   std::cout<<" Water Events = "<<water<<std::endl;
   std::cout<<" Carbon Events = "<<carbon<<std::endl;
   std::cout<<" Iron Events = "<<iron<<std::endl;
   std::cout<<" Lead Events = "<<lead<<std::endl;
   std::cout<<" Scintillator Events = "<<scintillator<<std::endl;
   std::cout<<"**********************************"<<std::endl;
   //	return variables;
}
//============================================================================================================================
// Main
