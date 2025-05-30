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
//using namespace globalV;
using namespace NUKECC_ANA;
//======================================================================
typedef VarLoop::Variable Var;
typedef Var2DLoop::Variable2D Var2D;
 
std::map<std::string, std::vector<CVUniverse*> > GetErrorBands(
    PlotUtils::ChainWrapper* chain) {
  typedef std::map<std::string, std::vector<CVUniverse*> > SystMap;

  SystMap error_bands;

  // CV
  error_bands[std::string("CV")].push_back(new CVUniverse(chain));


         if(RunCodeWithSystematics){
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
//
//void FillVariable(PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,bool isMC, int targetID, int targetZ, const string playlist, bool doDIS, std::vector<Var*>& variables,std::vector<Var2D*>& variables2d){

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

   Var* enu = new Var("Enu", "E#nu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrue);
   Var* ehad = new Var("Ehad", "Ehad (GeV)", Ehadbin, &CVUniverse::GetEhadGeV, &CVUniverse::GetEhadTrue);
   Var* Q2 = new Var("Q2", "Q2 (GeV^2)", Q2bin, &CVUniverse::GetQ2RecoGeV, &CVUniverse::GetQ2TrueGeV);
   Var* W = new Var("W", "W (GeV)", Wbin, &CVUniverse::GetWRecoGeV, &CVUniverse::GetWTrueGeV);
   Var* emu = new Var("Emu", "Emu (GeV)", Emubin, &CVUniverse::GetMuonEGeV, &CVUniverse::GetMuonETrueGeV);
   Var* x = new Var("x", "x", xbin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
   Var* y = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);
   Var* planeDNN = new Var("planeDNN", "planeDNN", planeDNNbin, &CVUniverse::GetplaneDNNReco, &CVUniverse::GetplaneDNNTrue);


   variables = {planeDNN,ehad,enu,emu, Q2, W, x, y}; 
   //variables = {planeDNN}; 
   
   Var2D* Q2_W = new Var2D(*Q2, *W);
   Var2D* enu_ehad = new Var2D(*enu, *ehad);
   Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
   variables2d = {Q2_W};//{enu_ehad, Q2_W};
   
   for (auto v : variables2d) v->InitializeAllHistograms(error_bands);
   for (auto v : variables) v->InitializeAllHistograms(error_bands,targetZ);

   int reco0=0;
   int reco00=0;
   int reco000=0;
   int reco1=0;
   int reco2=0;
   int reco3=0;
   int reco4=0; 
   int allcuts=0;

//------------------------------------------------Physics Sb stuff   
    sideBand transitionSB;
    sideBand continumSB;
    sideBand SB1;
    sideBand SB2;

    transitionSB.minQ2 = 1.0;
    transitionSB.maxQ2 = 200.0;
    transitionSB.minW = 1.5;
//    transitionSB.maxW = 2.0;
    transitionSB.maxW = 1.9;

    continumSB.minQ2 = 0.3;
//    continumSB.maxQ2 = 1.0;
    continumSB.maxQ2 = 0.9;
    continumSB.minW = 2.0;
    continumSB.maxW = 20.0;

    SB1.minQ2 = 1.0;
    SB1.maxQ2 = 200.0;
    SB1.minW = 1.9;
    SB1.maxW = 2.0;

    SB2.minQ2 = 0.9;
    SB2.maxQ2 = 1.0;
    SB2.minW = 2.0;
    SB2.maxW = 20.0;

//------------------------------------------------Physics Sb stuff   

   //=========================================
   // Targets combining Loop
   //=========================================
   for(int t = 0; t < targets.size(); t++){	  		
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
       std::vector<CVUniverse*> error_band_universes = band.second;
       for (auto universe : error_band_universes){
	 // Tell the Event which entry in the TChain it's looking at
	 universe->SetEntry(i);
	 //=========================================
	 // CUTS in each universe
	 //=========================================
	 if(!cutter->PassReco(universe,helicity)) continue;
	 reco1++;
	 if(!cutter->IsInMaterial(universe,targets[t],targetZ, /*anyTrakerMod*/false)) continue;
	 reco2++;
	 //if(!cutter->PassMuEnergyCut(universe)) continue;//We need these to be in the variable loop since you do not cut on Emu if you are plotting Emu -- ANF 2020.4.6
	 //if(!cutter->PassThetaCut(universe))continue;//Same thing for if we are plotting theta_mu. --ANF 2020.4.6

	 if(targets[t]<10 && universe->GetInt("NukeCC_targetID") != targets[t]) continue;
	 if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;	   
	 reco3++;
	 // if(universe->ShortName()=="cv") allcuts++;
         //if(doDIS){
         //  if(!cutter->PassDISCut(universe)) continue; 
         //}   

	 if(isMC){ 

	   for (auto v : variables2d){
	     if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(universe)) continue;
	     if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(universe)) continue;

             double recoQ2 = universe->GetQ2RecoGeV();
             double recoW = universe->GetWRecoGeV();
             double trueQ2 = universe->GetQ2TrueGeV();
             double trueW = universe->GetWTrueGeV();

             const bool isInTrans  = transitionSB.isInSideband(trueQ2, trueW);
             const bool isInContin = continumSB.isInSideband(trueQ2, trueW);
             const bool isInSB1 = SB1.isInSideband(trueQ2, trueW);
             const bool isInSB2 = SB2.isInSideband(trueQ2, trueW);
             const bool isTrueDIS  = cutter->PassTrueDISCut( universe );
           
//------------------------------------------------------Physics 2D scatter plots------------------------------------
             if(cutter->IsInTrueMaterial(universe,targets[t],targetZ,/*anyTrakerMod*/false)){
               if( cutter->passTrueCCQE(universe)){
                    v->m_selected_mc_q2.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}
	       else if( cutter->PassLowWCutTrue(universe)){
                    v->m_selected_mc_trueLowW.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}
	       else if( cutter->PassLowQ2TransTrue(universe)){
                    v->m_selected_mc_trueQ2.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}
	       else if( cutter->PasstransTrue(universe)){
                    v->m_selected_mc_truetrans.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}
	       else if( cutter->PassTrueContin( universe )){
                    reco0++;
                    v->m_selected_mc_lowQ2.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}
	       else if( cutter->PassTrueTrans( universe )){
                    reco00++;
                    v->m_selected_mc_lowW.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}
	       else if( cutter->PassTrueDISCut( universe )){
                    reco000++;
                    v->m_selected_mc_dis.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}
	     }
//------------------------------------------------------Physics 2D scatter plots------------------------------------
	   }
	   for (auto v : variables){
	     if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(universe)) continue;
	     if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(universe))continue;

             double recoQ2 = universe->GetQ2RecoGeV();
             double recoW = universe->GetWRecoGeV();
             double trueQ2 = universe->GetQ2TrueGeV();
             double trueW = universe->GetWTrueGeV();

             const bool isInTrans  = transitionSB.isInSideband(recoQ2, recoW);
             const bool isInContin = continumSB.isInSideband(recoQ2, recoW);
             const bool isTrueDIS  = cutter->PassTrueDISCut( universe );
             const bool isTrueinTrans  = (trueW < 2.0 ); 
             const bool isTrueinContin = (trueW > 2.0 && trueQ2 < 1.0);

//-------------------------------------------------------------Physics Sideband stuff---------------------
                if(isInTrans){

                    if(isTrueinTrans && ((cutter->passTrueCCQE(universe)) || (cutter->PassLowWCutTrue(universe)) || (cutter->PassLowQ2TransTrue(universe)) || (cutter->PasstransTrue(universe)))  ){
                      v->m_selected_mc_trans.GetComponentHist("True_Trans")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                   }
                    if(isTrueinContin && (cutter->PassLowQ2CutTrue(universe))){
                      v->m_selected_mc_trans.GetComponentHist("True_Contin")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                   }
                   if(isTrueDIS && (cutter->PassTrueDISCut(universe)) ){    
                      v->m_selected_mc_trans.GetComponentHist("True_DIS")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                   }
                }
                if(isInContin){
                    if(isTrueinTrans && ((cutter->passTrueCCQE(universe)) || (cutter->PassLowWCutTrue(universe)) || (cutter->PassLowQ2TransTrue(universe)) || (cutter->PasstransTrue(universe)))  ){
                      v->m_selected_mc_contin.GetComponentHist("True_Trans")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                   }
                    if(isTrueinContin && (cutter->PassLowQ2CutTrue(universe))){
                      v->m_selected_mc_contin.GetComponentHist("True_Contin")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                   }
                   if(isTrueDIS && (cutter->PassTrueDISCut(universe)) ){    
                      v->m_selected_mc_contin.GetComponentHist("True_DIS")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                   }

                }
                if(isTrueDIS){
                    v->m_selected_mc_sig.GetComponentHist("Signal")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                }
               if( !(cutter->IsInTrueMaterial(universe,targets[t],targetZ,/*anyTrakerMod*/false)) && targets[t]<10){
                 if( (universe->GetTargetZStart(targets[t]) > universe->GetVecElem("mc_vtx",2)) && universe->Var("planeDNN" ) > universe->GetTargetPlane(targets[t])-1 && universe->Var("planeDNN" ) < universe->GetTargetPlane(targets[t])+2 ){
                        if( isInTrans) {
                            v->m_hists_trans_mc_CH_US.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                        }
                        if( isInContin) {                                                                        
                            v->m_hists_contin_mc_CH_US.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                        }
                }

                 if( (universe->GetVecElem("mc_vtx",2) > universe->GetTargetZEnd(targets[t])) && universe->Var("planeDNN" ) > universe->GetTargetPlane(targets[t])-1 && universe->Var("planeDNN" ) < universe->GetTargetPlane(targets[t])+2){
                        if( isInTrans) {
                            v->m_hists_trans_mc_CH_DS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                        }
                        if( isInContin) {                                                                        
                            v->m_hists_contin_mc_CH_DS.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
                        }
                }
            }
//-------------------------------------------------------------Physics Sideband stuff---------------------
	   }

	 } else if(!isMC && !RunCodeWithSystematics){
	   for (auto v : variables2d){
	     if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(universe)) continue;
	     if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(universe)) continue;	     
	     v->m_selected_data_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe));
	   }
	   for (auto v : variables){
                double recoQ2 = universe->GetQ2RecoGeV();
                double recoW = universe->GetWRecoGeV();

                const bool isInTrans  = transitionSB.isInSideband(recoQ2, recoW);
                const bool isInContin = continumSB.isInSideband(recoQ2, recoW);

	     if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(universe)) continue;
	     if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(universe))continue;
             if(isInTrans){
                v->m_selected_data_reco_trans.univHist(universe)->Fill(v->GetRecoValue(*universe, 0));
             }
             if(isInContin){
                v->m_selected_data_reco_contin.univHist(universe)->Fill(v->GetRecoValue(*universe, 0));
             }

	   }
	 }
       } // End band's universe loop
     }// End Band loop
   }//End entries loop
 }//loop over targets

   std::cout<<"**********************************"<<std::endl;
   std::cout<<"Printing the ";
      isMC? std::cout<<"MC": std::cout<<"Data";
      std::cout<<" Summary "<<std::endl;
      std::cout<<" No cuts 1= "<<reco0<<std::endl;
      std::cout<<" No cuts 2= "<<reco00<<std::endl;
      std::cout<<" No cuts 3= "<<reco000<<std::endl;
   std::cout<<" Reco Cut = "<<reco1<<std::endl;
   std::cout<<" Material Cut = "<<reco2<<std::endl;
   std::cout<<" Plane/targetID Cuts = "<<reco3<<std::endl;
   std::cout<<" Muon Kinematics Cuts = "<<reco4<<std::endl;
   //   std::cout<<" Final CV= "<<allcuts<<std::endl; 
   std::cout<<"**********************************"<<std::endl;
   //	return variables;
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

cout<<__LINE__<<endl;	
   TString dir(argv[1]);
   int targetID = atoi(argv[2]);
   int targetZ = atoi(argv[3]);
   const string playlist= argv[4];
   
   bool doDIS=false;

cout<<__LINE__<<endl;	
   ROOT::Cintex::Cintex::Enable();
   TH1::AddDirectory(false);

   PlotUtils::DefaultCVUniverse::SetNFluxUniverses(100);
   PlotUtils::DefaultCVUniverse::SetNuEConstraint(false);
   PlotUtils::DefaultCVUniverse::SetAnalysisNuPDG(14);
   PlotUtils::DefaultCVUniverse::SetNonResPiReweight(true);

     
cout<<__LINE__<<endl;	
   NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(playlist);
   NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
   NukeCC_Binning  *binsDef = new NukeCC_Binning();
  
   PlotUtils::ChainWrapper* chainData = utils->GetChainWrapperDataPointer(playlist);
   PlotUtils::ChainWrapper* chainMC = utils->GetChainWrapperMCPointer(playlist);
   HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(playlist);
   
   double total_pot_data,total_pot_mc;
   utils->getPOT(total_pot_data,total_pot_mc);  
   double  MCscale=total_pot_data/total_pot_mc;
 
cout<<__LINE__<<endl;	
   std::vector<Var*> variablesMC,variablesData; 
   std::vector<Var2D*> variables2DMC,variables2DData; 

   TString histFileName;
   if(doDIS){
       histFileName = utils->GetHistFileName( "EventSelectionPhysicsDIS", FileType::kAny, targetID, targetZ, helicity ); 
     }
   else{
       histFileName = utils->GetHistFileName( "EventSelectionPhysics", FileType::kAny, targetID, targetZ, helicity );
     } 
 
cout<<__LINE__<<endl;	
  
   TFile fout(dir.Append(histFileName),"RECREATE");	
   
   utils->writePOT(fout);

cout<<__LINE__<<endl;	
   // For 1D variables 
   FillVariable(chainMC, helicity, utils, cutter,binsDef,true, targetID, targetZ, playlist,doDIS,variablesMC,variables2DMC);
       
   for (auto v : variablesMC) v->m_selected_mc_reco.SyncCVHistos();
   //for (auto v : variables2DMC) v->m_selected_mc_reco.SyncCVHistos();
     for (auto v : variables2DMC){
           v->m_selected_mc_reco.SyncCVHistos();
           v->m_selected_mc_lowW.SyncCVHistos();
           v->m_selected_mc_lowQ2.SyncCVHistos();
           v->m_selected_mc_dis.SyncCVHistos();
           v->m_selected_mc_q2.SyncCVHistos();
           v->m_selected_mc_trueLowW.SyncCVHistos();
           v->m_selected_mc_trueQ2.SyncCVHistos();
           v->m_selected_mc_truetrans.SyncCVHistos();
        }


cout<<__LINE__<<endl;	
   if(!RunCodeWithSystematics)FillVariable(chainData, helicity, utils, cutter,binsDef,false, targetID, targetZ, playlist,doDIS,variablesData,variables2DData);
   
cout<<__LINE__<<endl;	
   for (auto v : variablesData) v->m_selected_data_reco.SyncCVHistos();
//   for (auto v : variablesData) v->m_selected_data_reco_sb.SyncCVHistos();
   for (auto v : variablesData) v->m_selected_data_reco_trans.SyncCVHistos();
   for (auto v : variablesData) v->m_selected_data_reco_contin.SyncCVHistos();
   for (auto v : variables2DData) v->m_selected_data_reco.SyncCVHistos();
 
   for (auto v : variablesMC) {
     v->WriteAllHistogramsToFile(fout, true);
   }


cout<<__LINE__<<endl;	
   for (auto v : variablesData) {
     v->WriteAllHistogramsToFile(fout, false);
   }

   // Plotting If you want for 1D
      for(int i=0; i< variablesMC.size();i++){
      PlotStacked(variablesData[i]->m_selected_data_reco_trans.hist, variablesMC[i]->m_selected_mc_trans.GetHistArray(), MCscale, variablesMC[i]->m_selected_mc_trans.GetName(), variablesMC[i]->m_selected_mc_trans.GetName());
      PlotStacked(variablesData[i]->m_selected_data_reco_contin.hist, variablesMC[i]->m_selected_mc_contin.GetHistArray(), MCscale, variablesMC[i]->m_selected_mc_contin.GetName(), variablesMC[i]->m_selected_mc_contin.GetName());
      PlotStacked(variablesData[i]->m_selected_data_reco.hist, variablesMC[i]->m_selected_mc_sig.GetHistArray(), MCscale, variablesMC[i]->m_selected_mc_sig.GetName(), variablesMC[i]->m_selected_mc_sig.GetName());
cout<<__LINE__<<endl;	
   }//End 1D plotting 
   
   for (auto v : variables2DMC) {
     v->WriteAllHistogramsToFile(fout,true);
   }

   for (auto v : variables2DData) {
     v->WriteAllHistogramsToFile(fout,false);
}    
 
   //Plotting in 2D
   
   for(int i=0; i< variables2DMC.size();i++){
     Plot2D(variables2DMC[i]->m_selected_mc_reco.hist, variables2DMC[i]->GetName(), variables2DMC[i]->GetNameX(), variables2DMC[i]->GetNameY()); //Plotting line that I somehow cannot delete without producing memory errors, but no one else can reproduce. --ANF 2020.4.6
     //Plot2D(variables2DData[i]->m_selected_data_reco.hist, variables2DData[i]->GetName(), variables2DData[i]->GetNameX(),variables2DData[i]->GetNameY());
     
   }//End 2D plotting

}//End Main

