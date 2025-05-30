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
//#include "include/CVUniverse_faiza.h"
#include "include/CVUniverse.h"
#include "../../include/Variable_physicsSB.h"
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

void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID, int targetZ, const string playlist, bool doDIS=false);
 
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
// Main
int main(int argc, char *argv[]){
   ROOT::Cintex::Cintex::Enable();
   TH1::AddDirectory(false);
	
   if(argc==1){
     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
     std::cout<<"MACROS HELP:\n\n"<<
       "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n"<<
       "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
       "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" <<
       "\t-Playlist\t \t =\t eg. minervame1A"<< std::endl;
     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
     return 0;
   }

   TString dir(argv[1]);
   int targetID = atoi(argv[2]);
   int targetZ = atoi(argv[3]);
   const string playlist= argv[4];
   
   bool doDIS=false;
   const std::string mc_file_list("../../include/playlists/NukeCC_minervame1A_MC_Inextinguishable_merged.txt");
   const std::string data_file_list("../../include/playlists/NukeCC_minervame1A_DATA_Inextinguishable_merged.txt");
   const std::string plist_string(playlist);
   const std::string reco_tree_name("NukeCC");
   const bool wants_truth = false;
   const bool is_grid = false;

   PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth, is_grid);

    util.PrintMacroConfiguration("main");
  //=========================================
  // Systematics
  //=========================================

   PlotUtils::DefaultCVUniverse::SetNFluxUniverses(100);
   PlotUtils::DefaultCVUniverse::SetNuEConstraint(false);
   PlotUtils::DefaultCVUniverse::SetAnalysisNuPDG(14);
   PlotUtils::DefaultCVUniverse::SetNonResPiReweight(true);
     
   NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
   NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
   NukeCC_Binning  *binsDef = new NukeCC_Binning();
    
   PlotUtils::ChainWrapper* chainData = util.m_data;
   PlotUtils::ChainWrapper* chainMC   = util.m_mc;
   HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);
   double DataPot=  util.m_data_pot; 
   double MCPot=  util.m_mc_pot; 
   double  MCscale=DataPot/MCPot;
   std::cout<<" MCScale= "<<MCscale<<std::endl; 
 
   std::vector<Var*> variablesMC,variablesData; 
   std::vector<Var2D*> variables2DMC,variables2DData; 

   TString histFileName;
   if(RunCodeWithSystematics){
       histFileName = utils->GetHistFileName( "PhysicsBackgd_with_SYS_FullDet_q2WfromBranch_ME1A", FileType::kAny, targetID, targetZ, helicity ); 
     }
   else{
       histFileName = utils->GetHistFileName( "PhysicsBackgd_without_SYS", FileType::kAny, targetID, targetZ, helicity );
     } 
   TFile fout(dir.Append(histFileName),"RECREATE");	
   
   // For 1D variables 
   FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true,targetID, targetZ, plist_string,doDIS);
       
   for (auto v : variablesMC){
        v->m_hists_trans_in_trans.SyncCVHistos();
        v->m_hists_contin_in_trans.SyncCVHistos();
        v->m_hists_signal_in_trans.SyncCVHistos();
        v->m_hists_trans_in_contin.SyncCVHistos();
        v->m_hists_contin_in_contin.SyncCVHistos();
        v->m_hists_signal_in_contin.SyncCVHistos();
	v->m_hists_trans_mc_CH_US.SyncCVHistos();
        v->m_hists_contin_mc_CH_US.SyncCVHistos();
	v->m_hists_trans_mc_CH_DS.SyncCVHistos();
        v->m_hists_contin_mc_CH_DS.SyncCVHistos();
   }
     for (auto v : variables2DMC){
           v->m_selected_mc_dis.SyncCVHistos();
           v->m_selected_mc_truetrans.SyncCVHistos();
           v->m_selected_mc_trueQ2.SyncCVHistos();
           v->m_selected_mc_trueLowW.SyncCVHistos();
           v->m_selected_mc_q2.SyncCVHistos();
        }

   FillVariable(chainData, helicity, utils, cutter,binsDef,variablesData,variables2DData,false,targetID, targetZ, plist_string,doDIS);
 
   for (auto v : variablesData){
        v->m_selected_data_reco_trans.SyncCVHistos();
        v->m_selected_data_reco_contin.SyncCVHistos();
        //v->m_selected_data_reco.SyncCVHistos();
      }
   //for (auto v : variables2DData) v->m_selected_data_reco.SyncCVHistos();
 
   for (auto v : variablesMC) {
     v->WriteAllHistogramsToFile(fout, true);
   }

   for (auto v : variablesData) {
     v->WriteAllHistogramsToFile(fout, false);
   }

/*   
   // Plotting If you want for 1D
      for(int i=0; i< variablesMC.size();i++){
       
      PlotStacked(variablesData[i]->m_selected_data_reco_trans.hist, variablesMC[i]->m_selected_mc_trans.GetHistArray(), MCscale, variablesMC[i]->m_selected_mc_trans.GetName(), variablesMC[i]->m_selected_mc_trans.GetName());
      PlotStacked(variablesData[i]->m_selected_data_reco_contin.hist, variablesMC[i]->m_selected_mc_contin.GetHistArray(), MCscale, variablesMC[i]->m_selected_mc_contin.GetName(), variablesMC[i]->m_selected_mc_contin.GetName());
      PlotStacked(variablesData[i]->m_selected_data_reco.hist, variablesMC[i]->m_selected_mc_sig.GetHistArray(), MCscale, variablesMC[i]->m_selected_mc_sig.GetName(), variablesMC[i]->m_selected_mc_sig.GetName());

   }//End 1D plotting 
   for (auto v : variables2DMC) {
     v->WriteAllHistogramsToFile(fout,true);
   }
*/
cout<<"End of the macro"<<endl;
}//End Main

//-----------------------------------------Adding for Physics SB----------------------
        struct sideBand{
            double minQ2;
            double maxQ2;
            double minW;
            double maxW;

            bool isInSideband(double recoQ2, double recoW){
                return recoQ2 >= minQ2 && recoQ2 < maxQ2 && recoW >= minW && recoW < maxW;
            }
        };
//-----------------------------------------Adding for Physics SB----------------------

void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){

  std::map<std::string, std::vector<CVUniverse*> > error_bands =GetErrorBands(chain);

  std::vector<double> Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin;
  std::vector<double> planeDNNbin;

   if (doDIS){
     Enubin = binsDef->GetDISBins("Enu"); 
     Emubin = binsDef->GetDISBins("Emu"); 
     Ehadbin = binsDef->GetDISBins("Ehad");
     xbin = binsDef->GetDISBins("x");
     ybin = binsDef->GetDISBins("y");
     planeDNNbin = binsDef->GetDISBins("planeDNN");
     }
   else{
     Enubin = binsDef->GetEnergyBins("Enu"); 
     Emubin = binsDef->GetEnergyBins("Emu"); 
     Ehadbin = binsDef->GetEnergyBins("Ehad");
     xbin = binsDef->GetEnergyBins("x");
     ybin = binsDef->GetEnergyBins("y");
     planeDNNbin = binsDef->GetEnergyBins("planeDNN");
   }
   Q2bin = binsDef->GetSidebandBins("Q2");
   Wbin = binsDef->GetSidebandBins("W");

   Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
   Var* ehad = new Var("Ehad", "Ehad (GeV)", Ehadbin, &CVUniverse::GetEhadGeV, &CVUniverse::GetEhadTrueGeV);
   Var* Q2 = new Var("Q2", "Q2 (GeV^2)", Q2bin, &CVUniverse::GetQ2Reco, &CVUniverse::GetQ2True);
   Var* W = new Var("W", "W (GeV)", Wbin, &CVUniverse::GetWReco, &CVUniverse::GetWTrue);
   Var* emu = new Var("Emu", "Emu (GeV)", Emubin, &CVUniverse::GetMuonEGeV, &CVUniverse::GetMuonETrueGeV);
   Var* x = new Var("x", "x", xbin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
   Var* y = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);
   Var* planeDNN = new Var("planeDNN", "planeDNN", planeDNNbin, &CVUniverse::GetplaneDNNReco, &CVUniverse::GetplaneDNNTrue);
   
   variables = {emu}; 
   
   Var2D* Q2_W = new Var2D(*Q2, *W);
   Var2D* enu_ehad = new Var2D(*enu, *ehad);
   Var2D* emu_ehad = new Var2D(*emu, *ehad); 
   variables2d = {emu_ehad};
   
    std::vector<int> targetIDs;
    
    if( targetID > 10){
        targetIDs.push_back(14);
        targetIDs.push_back(24);
        targetIDs.push_back(34);
        targetIDs.push_back(44);
        targetIDs.push_back(54);
        targetIDs.push_back(64);
        targetIDs.push_back(74);
        targetIDs.push_back(84);
        targetIDs.push_back(94);
    }
    else{

        if(targetZ==26){
            targetIDs.push_back(1);
            targetIDs.push_back(2);
            targetIDs.push_back(3);
            targetIDs.push_back(5);
        }
        if(targetZ==82){
            targetIDs.push_back(1);
            targetIDs.push_back(2);
            targetIDs.push_back(3);
            targetIDs.push_back(4);
            targetIDs.push_back(5);
        }
        if(targetZ==6){
            targetIDs.push_back(3);
        }
    }

   for (auto v : variables2d) v->InitializeAllHistograms(error_bands);
   for (auto v : variables) v->InitializeAllHistograms(error_bands,targetID,targetZ);

   int reco0=0;
   int reco1=0;
   int reco2=0;
   int reco3=0;
   int reco4=0; 
   int reco5=0; 
   int reco6=0; 
   int allcuts=0;
   int t1=0;
   int t2=0;
   int t3=0;
   int t4=0;
   int c1=0;
   int d1=0;
   int dt=0;
   int data_dis=0;
   int c2=0;
   int c3=0;
   int c4=0;
   int s1=0;
   int ust=0;
   int dst=0;
   int usc=0;
   int dsc=0;
   int data1=0;
   int data2=0;
   int tDdis=0; 
   int tDtrans=0; 
   int tDcontin=0; 
 
   CVUniverse *dataverse = new CVUniverse(chain,0);
    
    sideBand transitionSB;
    sideBand continumSB;

    transitionSB.minQ2 = 1.0;
    transitionSB.maxQ2 = 200.0;
    transitionSB.minW = 1.5;
    transitionSB.maxW = 1.9;

    //continumSB.minQ2 = 0.3;
    continumSB.minQ2 = 0.5;
    continumSB.maxQ2 = 0.9;
    continumSB.minW = 2.0;
    continumSB.maxW = 20.0;
   //=========================================
   // Entry Loop
   //=========================================
   //for(int target = 0;target<targetIDs.size();target++){    //target loop
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
	 if(!cutter->PassReco(universe,helicity)) continue;
	 reco1++;
	 if(!cutter->IsInMaterial(universe,targetID,targetZ, /*anyTrakerMod*/false)) continue;
	 reco2++;
	 if(targetID<10 && universe->GetInt("NukeCC_targetID") != targetID) continue;
         reco3++;
	 if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;	   
	 reco4++;

         const double wgt = universe->GetWeight();

	   for (auto v : variables){
	     if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(universe)) continue;
	     if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(universe))continue;
	     reco5++;

             double recoQ2 = universe->GetQ2RecoGeV();
             double recoW = universe->GetWRecoGeV();
             double trueQ2 = universe->GetQ2TrueGeV();
             double trueW = universe->GetWTrueGeV();

             const bool isInTrans  = transitionSB.isInSideband(recoQ2, recoW);
             const bool isInContin = continumSB.isInSideband(recoQ2, recoW);

             const bool isTrueDIS  = cutter->PassTrueDISCut( universe );

             const bool isTrueinTrans  = (trueW < 2.0) && ((trueW < 1.5 ) || (trueQ2 < 1.0 && 1.5 <= trueW && trueW < 2.0) || (trueQ2 >= 1.0 && 1.5 <= trueW && trueW < 2.0) || (cutter->passTrueCCQE(universe)));
             const bool isTrueinContin = (trueW > 2.0 && trueQ2 < 1.0) && (trueW >= 2.0 && trueQ2 < 1.0);

            if(cutter->PassDISCut(universe)){d1++;}
            if(cutter->PassTrueDISCut(universe)){dt++;}
	     reco6++;
//-------------------------------------------------------------Physics Sideband stuff---------------------
                if(isInTrans){
                    t1++;
                   if(isTrueinContin){
                      v->m_selected_mc_trans.GetComponentHist("True_Contin")->Fill(v->GetRecoValue(*universe, 0), wgt);
                      v->m_hists_contin_in_trans.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                     t3++;
                   }
                    if(isTrueinTrans){
                      v->m_selected_mc_trans.GetComponentHist("True_Trans")->Fill(v->GetRecoValue(*universe, 0), wgt);
                      v->m_hists_trans_in_trans.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                     t2++;
                   }
                   if(isTrueDIS){    
                      v->m_selected_mc_trans.GetComponentHist("True_DIS")->Fill(v->GetRecoValue(*universe, 0), wgt);
                      v->m_hists_signal_in_trans.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                     t4++;
                   }

                }
               if(isInContin){
                  c1++;
                   if(isTrueinContin){
                      v->m_selected_mc_contin.GetComponentHist("True_Contin")->Fill(v->GetRecoValue(*universe, 0), wgt);
                      v->m_hists_contin_in_contin.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                    c3++;
                   }
                    if(isTrueinTrans){
                      v->m_selected_mc_contin.GetComponentHist("True_Trans")->Fill(v->GetRecoValue(*universe, 0), wgt);
                      v->m_hists_trans_in_contin.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                    c2++;
                   }
                   if(isTrueDIS){    
                      v->m_selected_mc_contin.GetComponentHist("True_DIS")->Fill(v->GetRecoValue(*universe, 0), wgt);
                      v->m_hists_signal_in_contin.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                    c4++;
                   }
                }
                if(isTrueDIS){
                  s1++;
                      v->m_selected_mc_sig.GetComponentHist("Signal")->Fill(v->GetRecoValue(*universe, 0), wgt);
                      v->m_hists_true_signal.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                }
               if( !(cutter->IsInTrueMaterial(universe,targetID,targetZ,false)) && targetID<10){
                 if( (universe->GetTargetZStart(targetID) > universe->GetVecElem("mc_vtx",2)) && universe->Var("planeDNN" ) > universe->GetTargetPlane(targetID)-1 && universe->Var("planeDNN" ) < universe->GetTargetPlane(targetID)+2 ){     
                        if( isInTrans) { ust++;
                            v->m_hists_trans_mc_CH_US.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                        }
                        if( isInContin) {      usc++;                                                                  
                            v->m_hists_contin_mc_CH_US.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                        }
                }

                 if( (universe->GetVecElem("mc_vtx",2) > universe->GetTargetZEnd(targetID)) && universe->Var("planeDNN" ) > universe->GetTargetPlane(targetID)-1 && universe->Var("planeDNN" ) < universe->GetTargetPlane(targetID)+2){    
                        if( isInTrans) {  dst++;
                            v->m_hists_trans_mc_CH_DS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                        }
                        if( isInContin) {  dsc++;                                                                     
                            v->m_hists_contin_mc_CH_DS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                        }
                }
            }


/*
               if( !(cutter->IsInTrueMaterial(universe,targetIDs[target],targetZ,false)) && targetIDs[target]<10){
                 if( (universe->GetTargetZStart(targetIDs[target]) > universe->GetVecElem("mc_vtx",2)) && universe->Var("planeDNN" ) > universe->GetTargetPlane(targetIDs[target])-1 && universe->Var("planeDNN" ) < universe->GetTargetPlane(targetIDs[target])+2){     
                 if( (universe->GetTargetZStart(targetIDs[target]) > universe->GetVecElem("mc_vtx",2)) && universe->Var("planeDNN" ) > universe->GetTargetPlane(targetIDs[target])-1 && universe->Var("planeDNN" ) < universe->GetTargetPlane(targetIDs[target])+2){     
                        if( isInTrans) {  ust++; 
                            v->m_hists_trans_mc_CH_US.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                        }
                        if( isInContin) {  usc++;                                                               
                            v->m_hists_contin_mc_CH_US.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                        }
                }

                 if( (universe->GetVecElem("mc_vtx",2) > universe->GetTargetZEnd(targetIDs[target])) && universe->Var("planeDNN" ) > universe->GetTargetPlane(targetIDs[target])-1 && universe->Var("planeDNN" ) < universe->GetTargetPlane(targetIDs[target])+2){    
                        if( isInTrans) { dst++; 
                            v->m_hists_trans_mc_CH_DS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                        }
                        if( isInContin) {     dsc++;                                                            
                            v->m_hists_contin_mc_CH_DS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), wgt);
                        }
                }
            }
*/
//-------------------------------------------------------------Physics Sideband stuff---------------------
         } //End 1D variables loop

	   for (auto v : variables2d){
	     if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(universe)) continue;
	     if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(universe)) continue;

             double recoQ2 = universe->GetQ2RecoGeV();
             double recoW = universe->GetWRecoGeV();
             double trueQ2 = universe->GetQ2TrueGeV();
             double trueW = universe->GetWTrueGeV();

//------------------------------------------------------Physics 2D scatter plots------------------------------------------------------------------
             if(cutter->IsInTrueMaterial(universe,targetID,targetZ,false)){
               if( cutter->passTrueCCQE(universe)){
                    v->m_selected_mc_q2.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), wgt);
                }
	       else if( cutter->PassLowWCutTrue(universe)){
                    v->m_selected_mc_trueLowW.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), wgt);
                }
	       else if( cutter->PassLowQ2TransTrue(universe)){
                    v->m_selected_mc_trueQ2.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), wgt);
                }
	       else if( cutter->PasstransTrue(universe)){
                    v->m_selected_mc_truetrans.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), wgt);
                }
	       else if( cutter->PassTrueDISCut( universe )){
                    v->m_selected_mc_dis.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), wgt);
                }
	     }
//------------------------------------------------------Physics 2D scatter plots---------------------------------------------------------------------
         } //End 2D variables loop
       } // End band's universe loop
     }// End Band loop
   }//End isMC loop

else{
        dataverse->SetEntry(i);
         reco0++;
	 if(!cutter->PassReco(dataverse,helicity)) continue;
	 reco1++;
	 if(!cutter->IsInMaterial(dataverse,targetID,targetZ, false)) continue;
	 reco2++;
	 if(targetID<10 && dataverse->GetInt("NukeCC_targetID") != targetID) continue;
         reco3++;
	 if( dataverse->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;	   
         reco4++;

	   for (auto v : variables){
	     if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(dataverse)) continue;
	     if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(dataverse))continue;
             reco5++;

                double recoQ2 = dataverse->GetQ2RecoGeV();
                double recoW = dataverse->GetWRecoGeV();

                const bool isintrans  = transitionSB.isInSideband(recoQ2, recoW);
                const bool isincontin = continumSB.isInSideband(recoQ2, recoW);

         if(cutter->PassDISCut(dataverse)){data_dis++;}
             if(isintrans){
                  data1++;
                v->m_selected_data_reco_trans.hist->Fill(v->GetRecoValue(*dataverse, 0));
             }
             if(isincontin){   data2++;
                v->m_selected_data_reco_contin.hist->Fill(v->GetRecoValue(*dataverse, 0));
             }
             reco6++;

           }//End 1D var loop       
         } //End else 
   }//End entries loop
// }//End target loop 
for(auto band : error_bands){
    std::vector<CVUniverse*> band_universes = band.second;
    for(unsigned int i_universe = 0; i_universe < band_universes.size(); ++i_universe)
      delete band_universes[i_universe];
  } 
    
      delete dataverse;
/*
   std::cout<<"**********************************"<<std::endl;
   std::cout<<"Printing the "; isMC? std::cout<<"MC": std::cout<<"Data";
   std::cout<<" Summary "              <<std::endl;
   std::cout<<" Before any cut = "     <<reco0<<std::endl;
   std::cout<<" Reco Cut = "           <<reco1<<std::endl;
   std::cout<<" Material Cut = "       <<reco2<<std::endl;
   std::cout<<" Plane/targetID Cuts = "<<reco3<<std::endl;
   std::cout<<" ML prob Cuts = "       <<reco4<<std::endl;
   std::cout<<"**********************************"<<std::endl;
   std::cout<<" MC trans = "           <<reco5<<std::endl;
   std::cout<<" T in T = "             <<reco6<<std::endl;
   std::cout<<"**********************************"<<std::endl;
*/
   std::cout<<"**********************************"<<std::endl;
   std::cout<<"Printing the ";isMC? std::cout<<"MC": std::cout <<"Data";
   std::cout<<" Summary "              <<std::endl;
   std::cout<<" No cuts = "            <<reco0<<std::endl;
   std::cout<<" Reco Cut = "           <<reco1<<std::endl;
   std::cout<<" Material Cut = "       <<reco2<<std::endl;
   std::cout<<" TargetID Cuts = "      <<reco3<<std::endl;
   std::cout<<" ML prob Cuts = "       <<reco4<<std::endl;
   std::cout<<" checking1 = "          <<reco5<<std::endl;
   std::cout<<" checking2 = "          <<reco6<<std::endl;
   std::cout<<"**********************************"<<std::endl;
   std::cout<<" MC trans = "           <<t1<<std::endl;
   std::cout<<" Data trans = "         <<data1<<std::endl;
   std::cout<<" trans in trans = "     <<t2<<std::endl;
   std::cout<<" contin in trans = "    <<t3<<std::endl;
   std::cout<<" sig in trans = "       <<t4<<std::endl;
   std::cout<<"**********************************"<<std::endl;
   std::cout<<" MC contin = "          <<c1<<std::endl;
   std::cout<<" Data contin = "        <<data2<<std::endl;
   std::cout<<" trans in contin = "    <<c2<<std::endl;
   std::cout<<" contin in contin = "   <<c3<<std::endl;
   std::cout<<" sig in contin = "      <<c4<<std::endl;
   std::cout<<"**********************************"<<std::endl;
   std::cout<<" MC True DIS = "        <<dt<<std::endl;
   std::cout<<" MC DIS = "             <<d1<<std::endl;
   std::cout<<" Data DIS = "           <<data_dis<<std::endl;
   std::cout<<" signal = "             <<s1<<std::endl;
   std::cout<<"**********************************"<<std::endl;
   std::cout<<" US in trans = "        <<ust<<std::endl;
   std::cout<<" DS in trans = "        <<dst<<std::endl;
   std::cout<<" US in contin = "       <<usc<<std::endl;
   std::cout<<" DS in contin = "       <<dsc<<std::endl;
   std::cout<<"**********************************"<<std::endl;

}
//============================================================================================================================

