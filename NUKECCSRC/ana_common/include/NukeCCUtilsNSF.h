#ifndef MNV_NUKECCUTILSNSF_h
#define MNV_NUKECCUTILSNSF_h 1


//#include "GlobalIncludes.h"
//#include "GlobalParameters.h" 
//#include "CCQENuCutsNSF.h"
#include "RooUnfold/RooUnfoldBayes.h"
#include "RooUnfold/RooUnfoldSvd.h"
#include "RooUnfold/RooUnfoldTUnfold.h"
#include "RooUnfold/RooUnfoldInvert.h"
#include "RooUnfold/RooUnfoldBinByBin.h"
#include "RooUnfold/RooUnfoldResponse.h"
#include "RooUnfold/RooUnfold.h"

#include "../include/NukeCC_Cuts.h"
//#include "CCQENuBinning.h"
//#include "include/ComputeUtils.h"
//#include "include/NeutronBlob.h"
//#include "include/NeutronBlob.h"
//#include "include/NeutronBlobBinning.h"
//#include "include/NeutronBlobCuts.h"
#include "CVUniverse.h"
#include "PlotUtils/MinervaUniverse.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvH3D.h"
#include "PlotUtils/AnaBinning.h"



using namespace PlotUtils;

// HadronReweight
namespace PlotUtils{
  class MnvHadronReweight;
  class FluxReweighter;
}

namespace NUKECC_ANA{

  class NukeCCUtilsNSF {

    typedef std::map<std::string,std::vector<double>*> StrVecMap;
    //! string-MnvND type for persistent map definition  
    //outside event loop
    typedef std::map<std::string, HistWrapper<CVUniverse>**> StrHist1DWrapperMap;
    typedef std::map<std::string, Hist2DWrapper<CVUniverse>**> StrHist2DWrapperMap;

private:
       //Michel Weights 
      std::vector<double> michel_weights; 

      //Target Mass Weights 
      std::vector<double> target_mass_weights; 

     //Proton Tracking Efficiency Weights 
      std::vector<double> track_eff_weights;  
      bool use_merged_files;


public:
    //! Default Constructor
    
    NukeCCUtilsNSF();
    
    NukeCCUtilsNSF(string playlist);

    //! Constructor with parameters
    NukeCCUtilsNSF( bool use_mergedFiles,string playlist );

    NukeCCUtilsNSF( bool use_mergedFiles, bool useFluxConstraint );

    //! Default Destructor
    ~NukeCCUtilsNSF();

    //! singleton gettor
    static NukeCCUtilsNSF& Get();

    public:
    // Some units
    StrHist1DWrapperMap cvhistos1D;
    StrHist2DWrapperMap cvhistos2D;
    
    //cutter
    NukeCC_Cuts *cutter;
    
    int incoming_pdg;

    //set the pot profile....
    double C_global_used_pot_mc;
    double C_global_tot_pot_mc;



    public:
      void bookHistos( TFile* file, TH2D** h, string var_name );
    
    std::string SetPlaylist(std::string playlist);

   void getPOT(double& total_pot_data,double& total_pot_mc );
   void writePOT(TFile *f, double DataPOT, double MCPOT);

HelicityType::t_HelicityType GetHelicityFromPlaylist( std::string playlist ) const;


std::vector<std::string> GetStdPlaylists( HelicityType::t_HelicityType helicity ) const;

      TString GetXSecHistName( int targetID, int targetZ, std::string var, int nucleon = 0, std::string cutName = "", bool fine = false );
      //! Extract the GENIE XSec histogram from this file
      MnvH1D *GetXSecHist( TFile *f, int targetZ, const std::string& var, HelicityType::t_HelicityType helicity, int targetID = 0, const std::string& cutName = "", bool isoCorrect = false, bool fine = false );
      //! Get necessary histograms and make xsec histo for tracker normalization
      MnvH1D *GetXSecRecoHistNorm( const std::string& var, FileType::t_FileType type, MatPair mat, HelicityType::t_HelicityType helicity, bool isoCorrect, TFile *f1 );
      //! Get necessary histograms and make avg of xsec histos for tracker normalization
      MnvH1D *GetXSecRecoHistNormAvg( const std::string& var, FileType::t_FileType type, int i_targetZ, HelicityType::t_HelicityType helicity, bool isoCorrect, TFile *f1  );

      //! Get necessary histograms and make xsec histo
      MnvH1D *GetXSecRecoHist( const std::string& var, FileType::t_FileType type, MatPair mat, HelicityType::t_HelicityType helicity, bool isoCorrect, TFile *f1 );
      //! Get necessary histograms and make avg xsec histo for these materials
      MnvH1D *GetXSecRecoHistAvg( const std::string var, FileType::t_FileType type, MatPairs mats,  HelicityType::t_HelicityType helicity, bool isoCorrect, TFile *f1 );
      //! Get necessary histograms and make avg xsec histo for materials with his Z
      MnvH1D *GetXSecRecoHistAvg( const std::string& var, FileType::t_FileType type, int i_targetZ, HelicityType::t_HelicityType helicity, bool isoCorrect, TFile *f1  );
      //! What Material pairs do we normally look at?
      MatPairs GetStdMatPairs( int nFaux = 0 ) const;
      //! What Material pairs do we use for normalization?
      MatPairs GetStdTrackerNormPairs( int targetZ = 0 ) const;
      //! What Material pairs do we use for normalization?
      MatPairs GetStdTrackerNormPairs( MatPair mat ) const;
      //! What Material pairs do we use for normalization?
      MatPairs GetStdTrackerNormPairs( int targetID, int targetZ ) const;
      //! What Material pairs do we group together for special purposes
      MatPairs GetSpecialMatPairs( int targetZ )const;
      //! Get the isoscalar correction that removes differences due to neutron excess
      PlotUtils::MnvH1D * GetIsoscalarCorrection( int targetZ, std::string var, int targetID = 0, const std::string& cutString = "", bool fine = false ) const;
      MnvH1D* EfficiencyCorrect( MnvH1D *h, MnvH1D *eff ) const;
      //! Name of the file that contains free nucleon histograms
      TString GetFreeNucleonXSecFileName() const;
      //! Add up the freen and freep cross sections
      MnvH1D * GetFreeNucSumXSec( int targetZ, std::string var, int targetID = 0, std::string cutName = "", bool fine = false ) const;
      //! Add up the freen and freep cross sections and average folding over all targets
      MnvH1D * GetFreeNucSumXSecAvg( int targetZ, const std::string& var, const std::string& cutName = "", bool fine = false ) const;
      //! Multiply average n+p xsec by A
      MnvH1D * GetIsoXSec( int targetZ, std::string var, int targetID = 0, std::string cutName = "", bool fine = false) const;
      //! Multiply average n+p xsec by A, average over folded targets is desired
      MnvH1D * GetIsoXSecAvg( int targetZ, const std::string& var, const std::string& cutName = "", bool fine = false) const;
      //! What Z and N do I use to describe the nucleus/material with this Zin?
      bool GetNucleusZN( int Zin, double &Z, double &N ) const;
      
      //bool DivideByFlux( MnvH1D* h, MnvH1D* &hFluxDraw, std::string playlist, bool reweightFlux = true ) const;// , int targetZ = 0) const;
      bool DivideByFlux( MnvH1D* h, std::string playlist, bool reweightFlux = true ) const;// , int targetZ = 0) const;
      bool DivideByIntegratedFlux( MnvH1D* h, std::string playlist, bool reweightFlux = true) const;// , int targetZ = 0) const;
      double GetTotalScatteringCenters( int targetZ, bool isMC ) const;

      std::string GetTotalXSecString( int targetID, int targetZ ) const;
      std::string GetTotalXSecUnits( ) const; 
      std::string GetSignalString() const;

      //! What latex string do I use to write this cross section
      std::string GetXSecString( const std::string& var, const MatPair& mat ) const;
      std::string GetXSecString( const std::string& var, const int targetID, const int targetZ ) const;



/*! @brief Fold a histogram using a migration matrix
        @param[out] hFolded Your MnvH1D which now contains the result of folding
        @param[in] hGenerated The histogram to fold, in generated kinematics
        @param[in] var Variable being unfolded
        @param[in] targetID targetID of the sample
        @param[in] targetZ targetZ of the sample
        @param[in] helicity neutrino or antineutrino
        @return true if everything worked
       */
      bool FoldHisto(
          MnvH1D* hFolded,
          const MnvH1D *hGenerated,
          const std::string& var,
          int targetID,
          int targetZ,
          HelicityType::t_HelicityType helicity,
          bool addSys = true,
	  bool isDIS  = true
          ) const;
 
      bool FoldHistoNormalized(
          MnvH1D* hFolded,
          const MnvH1D *hGenerated,
          const std::string& var,
          int targetID,
          int targetZ,
          HelicityType::t_HelicityType helicity,
          bool addSys = true
          ) const;

      TString GetHistFileName( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ, HelicityType::t_HelicityType helicity) const;
      TString GetHistFileNamePlaylist( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ, HelicityType::t_HelicityType helicity, std::string playlist) const;
TString HistDir( bool forceDisk  = false , bool ignorePlaylist = false, bool needsbluearc = false )const;

       TString GetVariationTag() const;
std::string GetHelicityString( HelicityType::t_HelicityType helicity ) const;

TString  GetHistName( const std::string& histType, FileType::t_FileType fType, const std::string& var, int targetID, int targetZ ) const;
TString  GetHistName( const std::string& histType, FileType::t_FileType fType, const std::string& var) const;
 
std::string GetFileTypeString( FileType::t_FileType fType ) const;
  
 //getting the Chainwrapper pointers....
  
  //POT stuff...

      //! Do I unfold this variable?
      bool UnfoldVar( const std::string& var ) const;
      //! How many iterations of unfolding?
      int UnfoldIterations( const std::string& var ) ;
      //! Force unfold on all
      void SetForceUnfoldAll( bool yn = true );
      //! Force do not unfold on all
      void SetForceNoUnfoldAll( bool yn = true );
      //! How many Unfolding Iterations?
      void SetUnfoldIter( int iter );

      //! Name of the file that contains XSec histograms
      TString GetXSecFileName() const;
 
      std::string fPlaylist;

  };
}
#endif
