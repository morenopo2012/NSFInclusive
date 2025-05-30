#ifndef Systematics_h
#define Systematics_h

//==============================================================================
// Get Several standard MINERvA systematics
//==============================================================================

#include "../../../NUKECCSRC/ana_common/include/CommonIncludes.h"
#include "../../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "../../../NUKECCSRC/ana_common/include/LateralSystematics.h"
#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/GeantHadronSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/TargetMassSystematics.h"
//#include "PlotUtils/MLVertexSystematics.h"
//#include "PlotUtils/MonaSystematics.h"
using namespace NUKECC_ANA;

typedef std::map<std::string, std::vector<CVUniverse*> > SystMap;

// Error Bands
std::map<std::string, std::vector<CVUniverse*> > GetErrorBands(PlotUtils::ChainWrapper* chain) {
  
  // return map
  SystMap error_bands;

  // CV 
  error_bands[std::string("CV")].push_back(new CVUniverse(chain));

  if(RunCodeWithSystematics) {
  
    //========================================================================
    // FLUX
    //========================================================================
    SystMap flux_systematics = 
    PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain,CVUniverse::GetNFluxUniverses());
  
    error_bands.insert(flux_systematics.begin(), flux_systematics.end());

    //========================================================================
    // GENIE
    //========================================================================
    SystMap genie_systematics = PlotUtils::GetGenieSystematicsMap<CVUniverse>(chain);
    // change that true to a switch on do_nonrespi_tune
    error_bands.insert(genie_systematics.begin(), genie_systematics.end());
    // includes standard Genie + Rvx1pi (final state pion normalization) + FaCCQE

    //========================================================================
    // MnvTunes
    //======================================================================== 
  
    // RPA
    SystMap RPA_systematics = PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
    error_bands.insert(RPA_systematics.begin(),RPA_systematics.end());

    // 2p2h
    SystMap _2p2h_systematics = PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
    error_bands.insert(_2p2h_systematics.begin(),_2p2h_systematics.end());

  
    //========================================================================
    //  Muons
    //========================================================================

    // Muon momentum
    // MINOS
    SystMap muonminosP_systematics = PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
    error_bands.insert(muonminosP_systematics.begin(), muonminosP_systematics.end());
  
    // Muon momentum
    // MINERvA
    SystMap muonminervaP_systematics = PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
    error_bands.insert(muonminervaP_systematics.begin(), muonminervaP_systematics.end());

    //Muon Momentum Resolution
    SystMap muonP_resolutions = PlotUtils::GetMuonResolutionSystematicsMap<CVUniverse>(chain);
    error_bands.insert(muonP_resolutions.begin(),muonP_resolutions.end());
  
    //MINOS Efficiency
    SystMap minos_efficiency = PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
    error_bands.insert(minos_efficiency.begin(),minos_efficiency.end());

    //========================================================================
    // Detector
    //========================================================================
  
    SystMap angle_systematics = PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);
    ///////////////////////NSFDefaults::beamThetaX_Err,NSFDefaults::beamThetaY_Err);
    error_bands.insert(angle_systematics.begin(), angle_systematics.end());


    // Target Mass Systematics
    SystMap targetMass_systematics = PlotUtils::GetTargetMassSystematicsMap<CVUniverse>(chain);
    error_bands.insert(targetMass_systematics.begin(), targetMass_systematics.end());

    //========================================================================
    // Particle Response Systematics (Detector)
    //========================================================================

    const bool use_ID = true;
    const bool use_OD = true;
    std::string name_tag = "allNonMuonClusters";
    const bool use_neutron = false;
    const bool use_new = true;
    const bool use_proton = true;
    //const bool use_nucl = false;
    //const bool use_tracker = false;
    //const bool use_ecal = false;
   // const bool use_hcal = false;
    // ID, OD, nucl, tracker, ecal, hcal
    SystMap response_systematics = PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain,use_ID,use_OD,name_tag, use_neutron, use_new, use_proton); //, use_ID, use_OD, use_nucl, use_tracker, use_ecal, use_hcal);

    error_bands.insert(response_systematics.begin(), response_systematics.end());

    // Old implementation with p3
    //const bool use_neutron = false;
    //const bool use_new = true;
    //const bool use_proton = true;
    // Particle response
    // -> 2nd argument: NEUTRON
    // -> 3rd argument: use_new_part_response
    // -> 4th argument: PROTON
    //SystMap response_systematics = PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain, use_neutron, use_new, use_proton);
    //error_bands.insert(response_systematics.begin(), response_systematics.end());


    //========================================================================
    // GEANT hadrons with MnvHadronReweight
    //========================================================================
    SystMap geant_hadron_systematics = PlotUtils::GetGeantHadronSystematicsMap<CVUniverse>(chain);
    error_bands.insert(geant_hadron_systematics.begin(), geant_hadron_systematics.end());


    //========================================================================
    // ML Vertex uncertainty
    //========================================================================
//    SystMap ml_vertex_systematics = PlotUtils::GetMLVertexSystematicsMap<CVUniverse>(chain);
//    error_bands.insert(ml_vertex_systematics.begin(), ml_vertex_systematics.end());

    //========================================================================
    // MoNA uncertainty
    //========================================================================
//    SystMap mona_systematics = PlotUtils::GetMonaSystematicMap<CVUniverse>(chain);
//    error_bands.insert(mona_systematics.begin(), mona_systematics.end());

  }

  return error_bands;
}
#endif  // Systematics_h
