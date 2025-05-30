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
#include "../../NUKECCSRC/ana_common/include/CommonIncludes.h"
#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "../include/Variable.h"
#include "../include/playlists/playlists.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "../../NUKECCSRC/ana_common/include/LateralSystematics.h"
#include "mec_common.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>  // For std::max_element
#include <iterator>   // For std::distance
//#include "Cintex/Cintex.h"
#include "../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Cuts.h"
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
#include "../include/plotting_functions.h"
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

    TString dir(argv[1]); //get the path and the output rootfile name
    int targetID = atoi(argv[2]);//atoi: string to integer. targetID: 1,2,3,4,5
    int targetZ = atoi(argv[3]);// targetZ: 26,81

    const string playlist= argv[4]; //which playlist is going to be working with
    bool doDIS=false;
    //bool isSYS=true;
    const std::string mc_file_list( get_mc_files(playlist, targetID) ); // Oscar, definition from ana/include/playlists/playlists.h
    const std::string data_file_list( get_data_files(playlist) );
    const std::string reco_tree_name("MasterAnaDev");
    const bool wants_truth = true;
    const bool is_grid = false;

    PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, playlist, wants_truth); //Oscar go to MacroUtil.cxx in plot Utils line 56?

    util.PrintMacroConfiguration("main");//Oscar, is printing the information see MacroUtil.cxx with this name function
    //=========================================
    // Systematics
    //=========================================
    //std::map<std::string, std::vector<CVUniverse*> > error_bands =
    //GetErrorBands(util.m_mc);

    //PlotUtils::MinervaUniverse::SetNuEConstraint(true);
    PlotUtils::MinervaUniverse::SetPlaylist(playlist);
    PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
    PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
    //PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
    PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);
    PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);
    //PlotUtils::MinervaUniverse::RPAMaterials(true);

    NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(playlist); // Oscar, go ../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h
    NukeCC_Cuts     *cutter  = new NukeCC_Cuts(); // go to NUKECCSRC/ana_common/src/NukeCC_Cuts.cxx
    NukeCC_Binning  *binsDef = new NukeCC_Binning();

    PlotUtils::ChainWrapper* chainData = util.m_data;
    PlotUtils::ChainWrapper* chainMC = util.m_mc;
    HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(playlist);
    double DataPot=  util.m_data_pot;
    double MCPot=  util.m_mc_pot;
    double  MCscale=DataPot/MCPot;

    std::cout<<" MCScale= "<<MCscale<<std::endl;
    std::vector<Var*> variablesMC,variablesData;
    std::vector<Var2D*> variables2DMC,variables2DData;

    TString histFileName = utils->GetHistFileNamePlaylist( "EventSelection", FileType::kAny, targetID, targetZ, helicity, playlist );

    //Works good for the grid submission
    //TString histFileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ );

    TFile fout(dir.Append(histFileName),"RECREATE");
    std::cout<<"helicity: "<<helicity<<std::endl;
    // For 1D variables
    FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true,targetID, targetZ, playlist,doDIS);
    for (auto v : variablesMC) v->m_selected_mc_reco.SyncCVHistos();
    for (auto v : variables2DMC) v->m_selected_mc_reco.SyncCVHistos();
    for (auto v : variables2DMC) v->m_selected_mc_reco_Lead.SyncCVHistos();
    for (auto v : variables2DMC) v->m_selected_mc_reco_Carbon.SyncCVHistos();
    for (auto v : variables2DMC) v->m_selected_mc_reco_Other.SyncCVHistos();
    //for (auto v : variablesMC) v->m_selected_mc_reco.SyncCVHistos();

    for (auto v : variablesMC) {
      v->WriteAllHistogramsToFile(fout, true);
    }

    for (auto v : variables2DMC) {
      v->WriteAllHistogramsToFile(fout,true);
    }

    // DATA
    std::cout << "Processing Data and filling histograms" << std::endl;
    FillVariable(chainData, helicity, utils, cutter,binsDef,variablesData,variables2DData,false,targetID, targetZ, playlist,doDIS);
    for (auto v : variables2DData) v->m_selected_data_reco.SyncCVHistos();

    for (auto v : variables2DData) {
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
  std::vector<double> ThetaMuBin,Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin, Eavbin, q3bin;
  std::vector<double> planeDNNbin, vtx_z;

  if (doDIS){
    planeDNNbin = binsDef->GetDISBins("planeDNN");
  }
  else{
    Enubin = binsDef->GetEnergyBins("Enu");
    Emubin = binsDef->GetEnergyBins("Emu");
    Ehadbin = binsDef->GetEnergyBins("Ehad");
    Q2bin = binsDef->GetEnergyBins("Q2");
    Wbin = binsDef->GetEnergyBins("W");
    xbin = binsDef->GetEnergyBins("x");
    ybin = binsDef->GetEnergyBins("y");
    Eavbin=binsDef->GetEnergyBins("Eavail");
    //Eavbin_Low=binsDef->GetEnergyBins("Eavail_Low");
    q3bin =binsDef->GetEnergyBins("q3");
    planeDNNbin = binsDef->GetDISBins("planeDNN");
    vtx_z = binsDef->GetEnergyBins("vtxz");
    //KEbin = binsDef->GetEnergyBins("KE_preFSI");
  }
   //Q2bin = binsDef->GetSidebandBins("Q2");
   //Wbin = binsDef->GetSidebandBins("W");


  Var* thetaMu = new Var("GetThetamuDeg", "GetThetamuDeg (Degree)", ThetaMuBin, &CVUniverse::GetThetamuDeg, &CVUniverse::GetThetamuTrueDeg);
  Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
  Var* ehad = new Var("Ehad", "Ehad (GeV)", Ehadbin, &CVUniverse::GetEhadGeV, &CVUniverse::GetEhadTrueGeV);
  Var* Q2 = new Var("Q2", "Q2 (GeV^2)", Q2bin, &CVUniverse::GetQ2RecoGeV, &CVUniverse::GetQ2TrueGeV);
  Var* W = new Var("W", "W (GeV)", Wbin, &CVUniverse::GetWRecoGeV, &CVUniverse::GetWTrueGeV);

  Var* emu = new Var("Emu", "Emu (GeV)", Emubin, &CVUniverse::GetMuonEGeV, &CVUniverse::GetMuonETrueGeV);
  Var* x = new Var("x", "x", xbin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
  Var* y = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);
  Var* planeDNN = new Var("planeDNN", "planeDNN", planeDNNbin, &CVUniverse::GetplaneDNNReco, &CVUniverse::GetplaneDNNTrue);
  Var* vtxz = new Var("vtxz", "vtxz", vtx_z, &CVUniverse::GetVertexZNew, &CVUniverse::GetVertexZTrueNew);

  //Var* Eavail = new Var("Eavailable", "Eavail (GeV)", Eavbin, &CVUniverse::GetECALHCALAvEnergy, &CVUniverse::GetEavailTrue); //GetTrackerECALAvEnergy
  Var* Eavail = new Var("Eavailable", "Eavail (GeV)", Eavbin, &CVUniverse::GetTrackerECALAvEnergy, &CVUniverse::GetEavailTrue);
  Var* q3 = new Var("q3", "q3 (GeV)", q3bin, &CVUniverse::Getq3Reco, &CVUniverse::Getq3Truth);
  //Var* KE_FSI=new Var("KE_FSI", "KE_FSI (GeV)", KEbin, &CVUniverse::GetKEFSITruth, &CVUniverse::GetKEFSITruth);

  //Var* E_avail_Low = new Var("Available Energy","Available Energy (GeV)",);

  //variables = {emu, ehad, enu, x, y, Q2, planeDNN};//{enu,ehad};
  //variables = {emu, ehad, enu, Eavail, q3};//{enu,ehad};
  variables = {emu, ehad, enu,planeDNN, vtxz};//{enu,ehad};
  Var2D* W_Q2 = new Var2D(*W, *Q2);
  //Var2D* enu_ehad = new Var2D(*enu, *ehad);
  Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
  Var2D* x_y = new Var2D(*x, *y);  // y var
  Var2D* x_Q2 = new Var2D(*x, *Q2);  // y var
  Var2D* Eav_q3=new Var2D(*Eavail, *q3);

  //variables2d = {emu_ehad, x_Q2, x_y, Eav_q3};//{enu_ehad, Q2_W};
  variables2d = {Eav_q3};//{enu_ehad, Q2_W};

  for (auto v : variables) v->InitializeAllHistograms(error_bands);
  for (auto v : variables2d) v->InitializeAllHistograms(error_bands);

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
    //targetIDs.push_back(targetID);
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

  double no_cuts = 0;
  int reco0=0;
  int reco1=0;
  int reco2=0;
  int reco3=0;
  int reco4=0;
  int reco5=0;
  int reco6=0;
  int reco7=0;
  int allcuts=0;
  int data0=0;
  int data1=0;
  int data2=0;
  int data3=0;
  int data4=0;
  int data5=0;

  CVUniverse *dataverse = new CVUniverse(chain,0);

  bool anyTrackerMod = false;
  //=========================================
  // Targets combining Loop
  //=========================================
  ////for(int t = 0; t < targetIDs.size(); t++){
  //=========================================
  // Entry Loop
  //=========================================
  TH1D* hist_fate_NC = new TH1D("hist_fate_NC", "KE - No Interaction (NS)", 50, 0, 2000);
  TH1D* hist_fate_ChEx = new TH1D("hist_fate_ChEx", "KE - No ChEx (ChEx)", 50, 0, 2000);
  TH1D* hist_fate_Elas = new TH1D("hist_fate_Elas", "KE - No Elas (Elas)", 50, 0, 2000);
  TH1D* hist_fate_Abs = new TH1D("hist_fate_Abs", "KE - No Abs (Abs)", 50, 0, 2000);
  TH1D* hist_fate_nn_50 = new TH1D("hist_fate_nn_50", "KE - Abs nn", 50, 0, 2000);
  TH1D* hist_fate_pn_51 = new TH1D("hist_fate_pn_51", "KE - Abs pn", 50, 0, 2000);
  TH1D* hist_fate_pp_52 = new TH1D("hist_fate_pp_52", "KE - Abs pp", 50, 0, 2000);
  TH1D* hist_fate_PionP = new TH1D("hist_fate_PionP", "KE - No PionP (PionP)", 50, 0, 2000);
  TH1D* hist_fate_MultN = new TH1D("hist_fate_MultN", "KE - No MultN (MN)", 50, 0, 2000);
  TH1D* hist_fate_KnockOut = new TH1D("hist_fate_KnockOut", "KE - KnockOut (KnockOut)", 50, 0, 2000);


  std::cout<<"# of entries = "<<chain->GetEntries()<<std::endl;
  int counter =0;
  int counter1=0;
  int counter2=0;

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
        //if(!cutter->PassReco(universe,helicity)) continue;
        reco1++;
	      //if(!cutter->IsInMaterial(universe,targetID,targetZ, anyTrackerMod)) continue;
	      //if(!cutter->IsInMaterial(universe,targetIDs[t],targetZ, /*anyTrakerMod*/false)) continue;
        double z;
        z = universe->GetVecElem((universe->GetAnaToolName() + "_vtx").c_str(),2);
        double minZ = 8614.65; // first plane module 85
        double maxZ = 8793.86; //first plane module 89
        //double emu = universe->GetEmu();
        //std::cout<<"emu: "<<emu<<endl;
	      reco2++;
	      //if(targetID<10 && universe->GetInt((universe->GetAnaToolName() + "_ANN_targetID").c_str()) != targetID) continue;
	      //if(targetID<10 && universe->GetInt((universe->GetAnaToolName() + "_targetID").c_str()) != targetID) continue;
	      //if(targetIDs[t]<10 && universe->GetInt("NukeCC_targetID") != targetIDs[t]) continue;
        reco3++;
	      //if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
	      
        //Regions::idx::t_RegionIdx region = cutter->GetRegion(universe, targetIDs[target], targetZ);
        //Regions::idx::t_RegionIdx region = cutter->GetRegion(universe, targetID, targetZ);

        //MEC cuts
        bool cuts = cutter->passesAll(universe);
        if(!cuts) continue;
        reco4++;
	      double ene = universe->GetDouble("MasterAnaDev_E");
        //if (ene < 2000) {continue;}
	      //This is to make a cut in the interaction model

	      int mc_type;
	      mc_type = universe->GetInt("mc_intType");
        //if(mc_type != 1) continue; //QE
        //if(mc_type != 2) continue; //Resonant


        // in which channel is this event?
        Channels::idx::t_ChannelIdx chanIdx = cutter->GetChannel(universe, helicity);
        int targetA = universe->GetInt("mc_targetA");
        //std::cout<<"Before function"<<std::endl;

        auto FateNumpreFSI = universe->GetFateNum();
        //std::cout<<"After function"<<std::endl;
        //std::cout<<"targetA: "<<targetA<<std::endl;

        std::vector<int> FateNumV = std::get<0>(FateNumpreFSI); //FateNumber
        std::vector<int> preFSIV = std::get<1>(FateNumpreFSI); //Particle ID
        std::vector<double> preFSIEn = std::get<2>(FateNumpreFSI); //Energy preFSI
        std::vector<double> FDEn = std::get<3>(FateNumpreFSI); //Energy First daughter

        //Get Scale Factor for Physiscs Background
        /*int sideBand = -1;
        if( chanIdx > 0 && chanIdx < 5){
        sideBand = 0;
        }
        else if (chanIdx == 5){
          sideBand = 1;
        }
        */
        //int sideBand = -1;
        //if( universe->GetWTrueGeV() < 2.0){
        //    sideBand = 0;
        //}
        //else if (universe->GetQ2TrueGeV() < 1.0){
        //    sideBand = 1;
        //}

        int flag_mass = 0;
        counter2++;
        for (auto v : variables2d){
	        //if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(universe)) continue;
          if( v->GetNameX()=="Emu")reco5++;
	        //if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(universe)) continue;
          if( v->GetNameX()=="Emu")reco6++;
	        //if( !cutter->PassDISCut( universe ))continue;
          reco7++;

       	  v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());

          if (targetA == 12){
            v->m_selected_mc_reco_Carbon.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); //GetTrueValueX GetRecoValueX
            double eav = universe->GetEavailTrue();
	          //double ene = universe->GetDouble("MasterAnaDev_E");
	          if (eav < 0.30 && ene < 2000){
              //std::cout<<" z target: "<<universe->GetInt("mc_targetZ") << std::endl;
		          //std::cout<<"Eav: "<< eav<< std::endl;
	            //std::cout << "Energy: " << ene << std::endl;	
              counter++;
	          }
            if (eav < 0.30){	
              counter1++;
	          }
       
            //Code-1=========================================================================================================================================================
            //This code is to double counting
            //std::cout<<"entry: "<<i<<std::endl;
            //std::cout<<" size of vector: "<< FateNumV.size() << std::endl;
            for(int i=0;i<FateNumV.size();i++){
              //std::cout<<"Resonant: "<<mc_type<<" with particle:"<< preFSIV[i]<< " Energy: "<<preFSIEn[i]<<" FateNum: "<<FateNumV[i] <<std::endl;
              double pionWeight = 1.0;
              //std::cout<<"Fate "<< FateNumV[i] <<std::endl;

              //Make sure is a pion
              //if (preFSIV[i] != 111 && preFSIV[i] != 211 && preFSIV[i] !=-211) continue;

              if(preFSIV[i] == 2212 && flag_mass == 0){
                //std::cout<<"QE proton "<< preFSIEn[i] <<std::endl;
                preFSIEn[i] = preFSIEn[i] - 938.72; //This is to get the KE only
                flag_mass = 1;
              }
              if(preFSIV[i] == 2112 && flag_mass == 0){
                //std::cout<<"QE neutron "<< preFSIEn[i] <<std::endl;
                preFSIEn[i] = preFSIEn[i] - 940.72;
                flag_mass = 1;
              }
              if(preFSIV[i] ==-211 && flag_mass == 0){
                //std::cout<<"QE C pion "<< preFSIEn[i] <<std::endl;
                preFSIEn[i] = preFSIEn[i] - 139.57;
                flag_mass = 1;
              }
              if(preFSIV[i]==211 && flag_mass == 0){
                //std::cout<<"QE C pion "<< preFSIEn[i] <<std::endl;
                preFSIEn[i] = preFSIEn[i] - 139.57;
                flag_mass = 1;
              }
              if(preFSIV[i] == 111 && flag_mass == 0){
                //std::cout<<"QE N pion "<< preFSIEn[i] <<std::endl;
                preFSIEn[i] = preFSIEn[i] - 135.57;
                flag_mass = 1;
              }

              //if(preFSIV[i] != 2212 && preFSIV[i] != 2112 && preFSIV[i] != 111 && preFSIV[i] != 211 && preFSIV[i] != -211 ){
              //    std::cout<<"ID: "<< preFSIV[i] << std::endl;
              //}

              //std::cout<<"Value: "<< preFSIEn[i] <<std::endl;
              //if(preFSIV[i] == 111 || preFSIV[i]==211 || preFSIV[i] ==-211 ){
              //    pionWeight = universe->getResPionFateTableDec2024(preFSIEn[i], FateNumV[i]);
              //  pionWeight = universe->getResonantPionWeightJan2024(preFSIEn[i], FateNumV[i]);
              //  pionWeight = universe->getCarbonResonantPionWeightJan2024(preFSIEn[i], FateNumV[i]);
              //pionWeight = 1.0;
              //if (preFSIEn[i] < 150 &&  (FateNumV[i] == 1 ||  FateNumV[i] == 3)){ std::cout << "Energy: " << preFSIEn[i] << " Fate: "<< FateNumV[i] << " weight: "<<pionWeight <<std::endl; }
              //pionWeight = 1.0;
              //}

              if(preFSIV[i] == 2112 || preFSIV[i]==2212 ){
                //pionWeight = universe->getQEFateTableDec2024(preFSIEn[i], FateNumV[i]);
                //pionWeight = universe->getQEWeightJan2024(preFSIEn[i], FateNumV[i]);
                //pionWeight = universe->getCarbonQEWeightJan2024(preFSIEn[i], FateNumV[i]);
		            pionWeight = 1.0;
		          }
              /*
              if( preFSIEn[i] < 41){
                std::cout<<" Energy: "<<preFSIEn[i]<<std::endl;
                std::cout<<" ID: "<<preFSIV[i]<<std::endl;
                std::cout<<" Fate: "<< FateNumV[i] << std::endl;
                std::cout<<" weight: "<<pionWeight<<std::endl;

                int mc_er_nPart = universe->GetInt("mc_er_nPart");
                std::cout<<" size of vector: "<< FateNumV.size() << std::endl;
                std::cout<<"i: "<<i<<" Value: "<< FateNumV[i] <<std::endl;
                std::cout<<" i * "<< "status *" << " er_FD *" <<" er_LD *" <<" mc_intTyp " << " er_ID *" <<" mc_er_E *" <<std::endl;
                for(int k=0;k<mc_er_nPart;k++){

                int status = universe->GetVecElemInt("mc_er_status",k);
                int id = universe->GetVecElemInt("mc_er_ID",k);
                double E = universe->GetVecElem("mc_er_E",k);
                int fd = universe->GetVecElemInt("mc_er_FD",k);
                int ld = universe->GetVecElemInt("mc_er_LD",k);
                int type = universe->GetInt("mc_intType");

                std::cout<<" "<<k <<"  * "<< status<<"  *  "    <<fd<<"   *   " <<ld<<"    *    " <<type <<  "     *    "      << id<<"  *  " <<E<<std::endl;
                }

              }*/

              std::pair<double,double> FateWeight= universe->getFateWeightSingle(FateNumV[i]);
              double PionFateWeight = FateWeight.first;
              double PionEnergy = FateWeight.second;

              double weight = universe->GetWeight();
              double Tweight= universe->GetTruthWeight();
		  
              double FinalWeight = weight*pionWeight;
              double TFinalWeight = Tweight*pionWeight;

              //if(PionFateWeight != 1){
              //std::cout<<"Weight: "<< weight <<" pion weight: "<<PionFateWeight <<" final weight: "<< weight*PionFateWeight<<" Fate: "<<FateNumV[i] << " preFSI: "<<preFSIV[i] <<" PionEnergy: "<<PionEnergy <<std::endl;
              if(FateNumV[i] == 1){
                v->m_selected_mc_reco_NoScattering.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),FinalWeight );
                v->m_selected_mc_truth_reco_NoScattering.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),TFinalWeight );
                hist_fate_NC->Fill(preFSIEn[i],TFinalWeight); //,TFinalWeight
                /*
                //std::cout<<"Weight: "<< weight <<" pion weight: "<<PionFateWeight <<" final weight: "<< weight*PionFateWeight<<" Fate: "<<FateNumV[i] << " preFSI: "<<preFSIV[i] <<" PionEnergy: "<<PionEnergy <<std::endl;
                std::cout<<" z target: "<<universe->GetInt("mc_targetZ") <<" counter: "<<counter <<" counter1: "<<counter1<<" counter2: "<<counter2 <<" Fate size:"<< FateNumV.size()<<" FateEntryvector: "<<i <<" Fate value:"<<FateNumV[i]<< " Part ID: "<<preFSIV[i] <<std::endl;

                int mc_er_nPart = universe->GetInt("mc_er_nPart");
                std::cout<<" size of vector: "<< FateNumV.size() << std::endl;
                std::cout<<"i: "<<i<<" Value: "<< FateNumV[i] <<std::endl;
                std::cout<<" i * "<< "status *" << " er_FD *" <<" er_LD *" <<" mc_intTyp " << " er_ID *" <<" mc_er_E *" <<std::endl;
                for(int k=0;k<mc_er_nPart;k++){

                  int status = universe->GetVecElemInt("mc_er_status",k);
                  int id = universe->GetVecElemInt("mc_er_ID",k);
                  double E = universe->GetVecElem("mc_er_E",k);
                  int fd = universe->GetVecElemInt("mc_er_FD",k);
                  int ld = universe->GetVecElemInt("mc_er_LD",k);
                  int type = universe->GetInt("mc_intType");

                  std::cout<<" "<<k <<"  * "<< status<<"  *  "    <<fd<<"   *   " <<ld<<"    *    " <<type <<  "     *    "      << id<<"  *  " <<E<<std::endl;
                }
                */
              }
              if(FateNumV[i] == 2){
                v->m_selected_mc_reco_ChargeExt.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), FinalWeight );
                v->m_selected_mc_truth_reco_ChargeExt.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),TFinalWeight );
                hist_fate_ChEx->Fill(preFSIEn[i],TFinalWeight);
              } //,TFinalWeight
              
              if(FateNumV[i] == 3){
                v->m_selected_mc_reco_Elasticity.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),FinalWeight  );
                v->m_selected_mc_truth_reco_Elasticity.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), TFinalWeight );
                hist_fate_Elas->Fill(preFSIEn[i],TFinalWeight);
                //if (preFSIEn[i] < 150 &&  (FateNumV[i] == 1 ||  FateNumV[i] == 3)){ std::cout << "Energy: " << preFSIEn[i] << " Fate: "<< FateNumV[i] << " weight: "<<pionWeight <<std::endl; }
				      }
              
              if(FateNumV[i] == 5){
                v->m_selected_mc_reco_Absorption.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),FinalWeight  );
                v->m_selected_mc_truth_reco_Absorption.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),TFinalWeight  );
                hist_fate_Abs->Fill(preFSIEn[i],TFinalWeight);}

              if(FateNumV[i] == 50){
                v->m_selected_mc_reco_fate_nn_50.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),FinalWeight  );
                v->m_selected_mc_truth_reco_fate_nn_50.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),TFinalWeight  );
                hist_fate_nn_50->Fill(preFSIEn[i],TFinalWeight);
              }
              
              if(FateNumV[i] == 51){
                v->m_selected_mc_reco_fate_pn_51.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),FinalWeight  );
                v->m_selected_mc_truth_reco_fate_pn_51.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),TFinalWeight  );
                hist_fate_pn_51->Fill(preFSIEn[i],TFinalWeight);
              }

              if(FateNumV[i] == 52){
                v->m_selected_mc_reco_fate_pp_52.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),FinalWeight  );
                v->m_selected_mc_truth_reco_fate_pp_52.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),TFinalWeight  );
                hist_fate_pp_52->Fill(preFSIEn[i],TFinalWeight);
              }

              if(FateNumV[i] == 8){
                v->m_selected_mc_reco_PionProd.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), FinalWeight );
                v->m_selected_mc_truth_reco_PionProd.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), TFinalWeight );
                hist_fate_PionP->Fill(preFSIEn[i],TFinalWeight);}
              
              if(FateNumV[i] == 9){
                v->m_selected_mc_reco_MultNuc.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), FinalWeight );
                v->m_selected_mc_truth_reco_MultNuc.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),TFinalWeight  );
                hist_fate_MultN->Fill(preFSIEn[i],TFinalWeight);
              }
              
              if(FateNumV[i] == 10){
                v->m_selected_mc_reco_KnockOut.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), FinalWeight );
                v->m_selected_mc_truth_reco_KnockOut.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), TFinalWeight );
                hist_fate_KnockOut->Fill(preFSIEn[i],TFinalWeight);


                //std::cout<<" size of vector: "<< FateNumV.size() << std::endl;
                //std::cout<<"i: "<<i<<" Value: "<< FateNumV[i] <<" fate: "<<preFSIV[i] <<std::endl;
                //std::cout<<"i: "<<i+1<<" Value: "<< FateNumV[i+1] <<" fate: "<<preFSIV[i+1] <<std::endl;
                /*
                int mc_er_nPart = universe->GetInt("mc_er_nPart");
                std::cout<<" size of vector: "<< FateNumV.size() << std::endl;
                std::cout<<"i: "<<i<<" Value: "<< FateNumV[i] <<std::endl;
                std::cout<<" i * "<< "status *" << " er_FD *" <<" er_LD *" <<" mc_intTyp " << " er_ID *" <<" mc_er_E *" <<std::endl;
                for(int k=0;k<mc_er_nPart;k++){

                  int status = universe->GetVecElemInt("mc_er_status",k);
                  int id = universe->GetVecElemInt("mc_er_ID",k);
                  double E = universe->GetVecElem("mc_er_E",k);
                  int fd = universe->GetVecElemInt("mc_er_FD",k);
                  int ld = universe->GetVecElemInt("mc_er_LD",k);
                  int type = universe->GetInt("mc_intType");

                  std::cout<<" "<<k <<"  * "<< status<<"  *  "    <<fd<<"   *   " <<ld<<"    *    " <<type <<  "     *    "      << id<<"  *  " <<E<<std::endl;
                }

                */

               }
            } //end of FateNumV.size()

            //Code-1=========================================================================================================================================================


            /*
            //Code-2=========================================================================================================================================================
            //This code is for the most energetic

            // Find the iterator pointing to the largest element
            auto maxElement = std::max_element(preFSIEn.begin(), preFSIEn.end());
            int index;

            // Check if the vector is not empty
            if (maxElement != preFSIEn.end()) {
              // Get the index of the largest element
               ndex = std::distance(preFSIEn.begin(), maxElement);
              //std::cout << "The largest value is: " << *maxElement << " at index: " << index << std::endl;
            } else {
              continue;
            }


            double weight = universe->GetWeight();
            double Tweight= universe->GetTruthWeight();

            double pionWeight = 1.0;
            if(preFSIV[index] == 111 || preFSIV[index]==211 || preFSIV[index] ==-211 ){
              pionWeight = universe->getFateTable(preFSIEn[index], FateNumV[index]);
            }

            double FinalWeight = weight*pionWeight;
            double TFinalWeight = Tweight*pionWeight;


            //if(PionFateWeight != 1){
            //std::cout<<"Weight: "<< weight <<" pion weight: "<<PionFateWeight <<" final weight: "<< weight*PionFateWeight<<" Fate: "<<FateNumV[i] << " preFSI: "<<preFSIV[i] <<" PionEnergy: "<<PionEnergy <<std::endl;
            if(FateNumV[index] == 1){
              v->m_selected_mc_reco_NoScattering.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),FinalWeight );
              v->m_selected_mc_truth_reco_NoScattering.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),TFinalWeight );
              //std::cout<<"Weight: "<< weight <<" pion weight: "<<PionFateWeight <<" final weight: "<< weight*PionFateWeight<<" Fate: "<<FateNumV[i] << " preFSI: "<<preFSIV[i] <<" PionEnergy: "<<PionEnergy <<std::endl;
            }
            if(FateNumV[index] == 2){
              v->m_selected_mc_reco_ChargeExt.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), FinalWeight );
              v->m_selected_mc_truth_reco_ChargeExt.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),TFinalWeight  );
            }
            
            if(FateNumV[index] == 3){
              v->m_selected_mc_reco_Elasticity.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),FinalWeight*0.5  );
              v->m_selected_mc_truth_reco_Elasticity.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), TFinalWeight*0.5 );
            }           
            if(FateNumV[index] == 5){
              v->m_selected_mc_reco_Absorption.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),FinalWeight  );
              v->m_selected_mc_truth_reco_Absorption.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),TFinalWeight  );
            }
            if(FateNumV[index] == 8){
              v->m_selected_mc_reco_PionProd.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), FinalWeight );
              v->m_selected_mc_truth_reco_PionProd.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), TFinalWeight );}
            if(FateNumV[index] == 9){
              v->m_selected_mc_reco_MultNuc.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), FinalWeight );
              v->m_selected_mc_truth_reco_MultNuc.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),TFinalWeight  );}
            if(FateNumV[index] == 10){
              v->m_selected_mc_reco_KnockOut.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), FinalWeight );
              v->m_selected_mc_truth_reco_KnockOut.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), TFinalWeight );}

             //Code-2=========================================================================================================================================================
            */




          }//end of "if(targetA == 12)" or 207

          if (targetA == 207){
            v->m_selected_mc_reco_Lead.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}
          if (targetA != 207 && targetA != 12 ){
            v->m_selected_mc_reco_Other.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}

	      }// end of "for (auto v : variables2d)"

	      for (auto v : variables){
	        //if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(universe)) continue;
	        if( v->GetName()=="Emu")reco6++;
          // if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(universe))continue;
	        if (v->GetName()=="Enu") reco7++;
	        //if( !cutter->PassDISCut( universe )) continue;
	        v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
          /*if (targetID < 10){
            if(cutter->IsInTrueMaterial(universe,targetID,targetZ, anyTrackerMod) && cutter->PassTrueFiducial(universe)){
              if (cutter->PassTrueDISCut( universe ) && cutter->PassTrueMuEnergyCut(universe) && cutter->PassTrueThetaCut(universe))
              v->m_selected_mc_reco_true_target.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
            }
            else{
              if (false ) //region ==3
              v->m_selected_mc_reco_USCH.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              else if (false) //region ==4
              v->m_selected_mc_reco_DSCH.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
              else
              v->m_selected_mc_reco_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());

            }//Close else condtion
          }*/

          //else {
            //if (!(cutter->IsInTrueMaterial(universe,targetID,targetZ, anyTrackerMod) && cutter->PassTrueFiducial(universe)))
            //v->m_selected_mc_reco_other.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
          //}//Close targetID<10 if condtion
          //if ( cutter->IsInTrueMaterial(universe,targetID,targetZ, anyTrackerMod)){
            // if (cutter->PassTrueDISCut( universe )) {
            // if(sideBand==0) v->m_selected_trans_SB.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
            // if(sideBand==1) v->m_selected_contin_SB.univHist(universe)->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
          //}
	        // v->m_selected_mc_sb.GetComponentHist("MC")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
	        //if( cutter->PassDISCut( universe ))
	        // v->m_selected_mc_sb.GetComponentHist("DIS")->Fill(v->GetRecoValue(*universe, 0), universe->GetWeight());
	      }
      } // End band's universe loop
     }// End Band loop
    }//End of isMC

    if(!isMC){
      dataverse->SetEntry(i);

      //MEC cuts
      bool cuts = cutter->passesAll(dataverse);
      if(!cuts) continue;

      data0++;
	    //if( dataverse->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
      data4++;
	    //if(!cutter->PassReco(dataverse,helicity)) continue;
	    data1++;
	    //if(!cutter->IsInMaterial(dataverse,targetID,targetZ, anyTrackerMod)) continue;
	    //if(!cutter->IsInMaterial(dataverse,targetIDs[t],targetZ, false)) continue;
      double z_data;
      z_data = dataverse->GetVecElem((dataverse->GetAnaToolName() + "_vtx").c_str(),2);
      double minZ = 8614.65; // first plane module 85
      double maxZ = 8793.86; //first plane module 89

	    data2++;
	    //if(targetID<10 && dataverse->GetInt((dataverse->GetAnaToolName() + "_ANN_targetID").c_str()) != targetID) continue;
	    //if(targetID<10 && dataverse->GetInt((dataverse->GetAnaToolName() + "_targetID").c_str()) != targetID) continue;
	    //if(targetIDs[t]<10 && dataverse->GetInt("NukeCC_targetID") != targetIDs[t]) continue;
      data3++;
	    //if( dataverse->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
      //data4++;
 	    for (auto v : variables2d){
	     //if( v->GetNameX()!="Emu" && v->GetNameY()!="Emu")  if(!cutter->PassMuEnergyCut(dataverse)) continue;
	     //if( v->GetNameX()!="ThetaMu" && v->GetNameY()!="ThetaMu")  if(!cutter->PassThetaCut(dataverse)) continue;
	     //if( !cutter->PassDISCut(dataverse)) continue;
	     v->m_selected_data_reco.hist->Fill(v->GetRecoValueX(*dataverse), v->GetRecoValueY(*dataverse));
	    }
	    for (auto v : variables){
	      //if( v->GetName()!="Emu")   if(!cutter->PassMuEnergyCut(dataverse)) continue;
	      //if( v->GetName()!="ThetaMu") if(!cutter->PassThetaCut(dataverse))continue;
	      //if (v->GetName()=="Enu") data5++;
	      //if( !cutter->PassDISCut(dataverse)) continue;
        data5++;
	      v->m_selected_data_reco.hist->Fill(v->GetRecoValue(*dataverse, 0));
	      //v->m_selected_data_reco_sb.hist->Fill(v->GetRecoValue(*dataverse, 0));
      }
    }

	  if(i%500000==0) std::cout << (i/100) << " hello 6k " << std::endl; //Eraseme
  }//End entries loop

 //}//Target loop closed

  std::cout<<" with energy lower than 2000:" << counter <<std::endl;
  std::cout<<" Total entries:" << counter1 <<std::endl;

  TFile *out_file = new TFile("KE_PB.root","RECREATE");
  hist_fate_NC->Write();
  hist_fate_ChEx->Write();
  hist_fate_Elas->Write();
  hist_fate_Abs->Write();
  hist_fate_nn_50->Write();
  hist_fate_pn_51->Write();
  hist_fate_pp_52->Write();
  hist_fate_PionP->Write();
  hist_fate_MultN->Write();
  hist_fate_KnockOut->Write();
  out_file->Close();

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
   std::cout<<" No cuts = "<<data0<<std::endl;
   std::cout<<" ML Cut = "<<data4<<std::endl;
   std::cout<<" Reco Cut = "<<data1<<std::endl;
   std::cout<<" Material Cut = "<<data2<<std::endl;
   std::cout<<" Plane/targetID Cuts = "<<data3<<std::endl;
   std::cout<<" DIS cut = "<<data5<<std::endl;
   std::cout<<"**********************************"<<std::endl;
   //	return variables;
}
//============================================================================================================================
// Main
