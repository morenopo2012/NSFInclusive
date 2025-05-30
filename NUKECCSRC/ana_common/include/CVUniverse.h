// ========================================================================
// Base class for an un-systematically shifted (i.e. CV) universe.  Implement
// your own base class with the functions you need. I've implemented GetEnu(),
// GetMuonE() and GetHadronE() as examples: you'll have other variables you
// want.
//
// To add a systematic, inherit from this class and override whichever
// functions you need to. For a "vertical" error, this will mean overriding the
// GetWeight() function to modify the event weight. For a "lateral" error, this
// will mean overriding the function that calculates the quantity that is being
// shifted (muon energy, or hadronic energy or whatever).
//
// For examples of each of those two cases, see ./LateralSystematics.h and
// PlotUtils/GenieSystematics.h. For an example of how to put the whole thing
// together and actually *use* the classes, see the runEventLoop.C macro in
// this directory. `root -l -b load.C+ runEventLoop.C+`
// ========================================================================
#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

//#include "NukeCCsetbranchaddrs.h"
#include "CommonIncludes.h"
#include "PlotUtils/MinervaUniverse.h"
#include "PlotUtils/ChainWrapper.h"
#include <iostream>
//#include "CCQENuUtilsNSF.h"
//using namespace globalV;

#include "PlotUtils/CaloCorrection.h"
//#include "NukeCCvars.h"
namespace NUKECC_ANA
{
class CVUniverse : public PlotUtils::MinervaUniverse {
  public:
    #include "PlotUtils/MuonFunctions.h" // GetMinosEfficiencyWeight
    #include "PlotUtils/TruthFunctions.h" //Getq3True
    #include "PlotUtils/WeightFunctions.h"
    #include "PlotUtils/RecoilEnergyFunctions.h"

    // Constructor
    CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma=0);
    CVUniverse();

    // Destructor
    virtual ~CVUniverse();
    // All functions we write here should be 'virtual', so that the universes
    // that inherit from CVUniverse can override them.

    // ========================================================================
    // Get Weight
    //
    // We override the various weight getting functions herein in different
    // vertical systematic universe classes.
    // ========================================================================


  virtual double GetWeight() const;
  virtual double GetWeightEmu() const;
  virtual double GetWeightQ2() const;
  virtual double GetWeighty() const;
  //virtual double GetLowQ2PiWeight(std::string channel, int shift) const;
  virtual double GetTruthWeight() const;
  virtual Long64_t GetMaxEntries();
  virtual double GetMuonE() const { return GetVecElem("NukeCC_leptonE",3); }

  //      virtual void  Show(Long64_t entry = -1);


  virtual double calcq3(const double Q2,const double Enu,const double Emu	)const;
  virtual double calcq0(const double Enu,const double Emu	)const;
  //! Turn off branches that I don't use
  //virtual bool SetBranches( bool useWeights = true );
      //! Get the Region idx used for the region vtx/planecode stacks plots
      Regions::idx::t_RegionIdx GetRegion( int targetID, int targetZ );

  virtual double GetMuonCurve()  const;
  virtual double GetFiducial()  const;
  virtual double GetTdead()  const;
  virtual double GetTargetID()  const;
  virtual double GetANNTargetID()   const;
  virtual double GetHelicity()  const;
  virtual double GetTrueHelicity()  const;

  virtual double GetECALHCALAvEnergy() const;
  virtual double calcECALHCALAvEn() const;

  virtual double GetTrackerECALAvEnergy() const;
  virtual double calcTrackerECALAvEn() const;
  //MECAna path to do it
  virtual double MECAnaGetTrackerECALAvEnergy() const;
  virtual double MECAnacalcTrackerECALAvEn() const;
  virtual double MECAnaGetq3Reco() const;
  virtual double MECAnacalcRecoq3() const;

  //Low Recoil
  virtual double GetEavalReco_Low() const;

  virtual std::tuple<std::vector<int>,std::vector<int>,std::vector<double>,std::vector<double> > GetFateNum() const;
  virtual int calFateNum() const;
  virtual std::tuple<std::vector<int>,std::vector<int>,std::vector<double>,std::vector<double> > calFateNumVector() const;
  virtual double getGenieBEinMeV(int A) const;
  virtual std::pair<std::vector<double>,std::vector<double> > getFateWeight(double FateNum) const;
  virtual double getFateTable(double piEnergy, int FateNum) const;
  virtual std::pair<double,double> getFateWeightSingle(double FateNum) const;

  virtual double getLeadResPionTable(double piEnergy, int FateNum) const;
  virtual double getLeadQETable(double piEnergy, int FateNum) const;
  virtual double getCarbonResPionTable(double piEnergy, int FateNum) const;
  virtual double getCarbonQETable(double piEnergy, int FateNum) const;

  virtual double getQEWeightJan2024( double piEnergy, int FateNum) const;
  virtual double getResonantPionWeightJan2024( double piEnergy, int FateNum) const;
  virtual double getCarbonQEWeightJan2024( double piEnergy, int FateNum) const;
  virtual double getCarbonResonantPionWeightJan2024( double piEnergy, int FateNum) const;

  //virtual double GetRecoilEnergy() const;
  virtual double GetEhadGeV() const;
  virtual double GetMuonEGeV() const;
  virtual double GetMuonETrueGeV() const;
  virtual double GetThetamuTrueDeg()  const ;
  //virtual double GetThetamu()  const ;
  virtual double GetMuonPt()  const ;
  virtual double GetMuonPz()  const ;
  virtual double GetlepPtTrue()  const ;
  virtual double GetQ2Reco()     const ;
  virtual double GetWReco()     const;
  virtual double GetxReco()     const;
  virtual double GetyReco()     const;
  virtual double GetplaneDNNReco() const;

  virtual double GetQ2RecoGeV()     const ;
  virtual double GetWRecoGeV()     const;
  virtual double GetEnu() const;
  virtual double GetEnuGeV() const;
  virtual double GetEnuTrueGeV() const;

  virtual double GetEavReco() const;
  virtual double GetEavailReco() const;
  virtual double GetEavailRecoECal() const;
  virtual double Getq3Reco() const;
  virtual double calcRecoq3() const;
  virtual double Getq0Truth() const;
  //virtual double GetEavTrue() const;
  virtual double GetEavailTrue() const;
  virtual double calcEavailTrue() const;
  virtual double Getq3Truth() const;
  virtual double calcTrueq3() const;
  virtual double Convert2PlaneID(double vtx) const;

  virtual double GetCVWeight() const;
  virtual double calcCVWeight() const;
  // for true variables
  virtual double GetTruthNuPDG() const;
  virtual double GetCurrent()    const;
  virtual double GetEhadTrue() const;
  virtual double GetEhadTrueGeV() const;
  virtual double GetThetamuTrue( )const;

  virtual double GetQ2IncTrue()     const ;
  virtual double GetWTrue()     const;
  virtual double GetyTrue()     const;
  virtual double GetxTrue()     const;
  virtual double GetplaneDNNTrue()   const;

  virtual double GetQ2TrueGeV()     const ;
  virtual double GetWTrueGeV()     const;
  virtual double GetThetamuDeg()  const;
  //virtual double Getq3True()const;
  //virtual double Getq0True() const;
  virtual double GetVertexZNew() const;
  virtual double GetVertexZTrueNew() const;
    double virtual GetVertexZMy() const;
    double virtual GetVertexZTrueMy() const;

  virtual double GetCalRecoilEnergy() const; // added for new particle response sys.
  virtual double GetNonCalRecoilEnergy() const; //added for new particle response sys.
  virtual double ApplyCaloTuning(double calRecoilE) const;
  
  virtual double GetMuonPZ() const;
  virtual double GetMuonPZTrue() const;

  virtual double calcQ2(const double Enu,const double Emu,const double Thetamu)const;
  virtual double calcW(const double Q2,const double Ehad) const;
  virtual double calcX(const double Q2,const double Enu, const double Emu)  const ;
  virtual double calcY(const double Enu,const double Ehad)const;


virtual double Var( const std::string& var, bool useTrue = false );

  virtual  int GetTargetFromSegment( int segment, int& vtx_module, int& vtx_plane ) const ;
  virtual double GetVtxEnergy(  unsigned int shift = VTX_BLOB_SHIFT );

  virtual double GetCCQERecoil( unsigned int shift = VTX_BLOB_SHIFT );

  virtual double GetMuonP() const;

      virtual int GetAtomicNumber();
      //! Map Atomic Number to Its Atomic Mass
      virtual int AtomicNumberToMass( int targetZ ) const;
      //! Map Target to the Minimum Plane Number Bin
      virtual int GetTargetMinBin( int targetID ) const;
      //! Map Target to the Maximum Plane Number Bin
      virtual int GetTargetMaxBin( int targetID ) const;
      //! Map Target to its Plane Number
      virtual int GetTargetPlane( int targetID ) const;
      //! Map Plane to its respective Target ID
      virtual int GetPlaneTargetID( int plane ) const;
      //! Map Target to the Plane immediately upstream of it
      virtual int GetTargetUSPlane( int targetID ) const;
      //! Map Target to the Plane immediately downstream of it
      virtual int GetTargetDSPlane( int targetID ) const;
      //! Map Target to the Plane immediately downstream of it
      virtual double GetTargetZStart( int targetID ) const;
      //! Map Target to the Plane immediately downstream of it
      virtual double GetTargetZEnd( int targetID ) const;

    // ARACHNE EVENT DISPLAY
    // data
    virtual int GetRunN() const;
    virtual int GetSubRunN() const;
    virtual int GetGateN() const;
    virtual int GetSliceN() const;
    virtual int Getvtx0N() const;
    virtual int Getvtx1N() const;
    virtual int Getvtx2N() const;
    virtual int Getvtx3N() const;

    // MC
    virtual int GetMCRunN() const;
    virtual int GetMCSubRunN() const;
    virtual int GetMCGateN() const;
    virtual int GetMCSliceN() const;
};
}

#endif
