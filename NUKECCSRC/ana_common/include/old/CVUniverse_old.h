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

//#include "NukeCCvars.h"
//#include "NukeCCsetbranchaddrs.h"
#include "CommonIncludes.h"
#include "PlotUtils/DefaultCVUniverse.h"
#include "PlotUtils/ChainWrapper.h"
#include <iostream>
//#include "NukeCC_Cuts.h"


namespace NUKECC_ANA
{
class CVUniverse : public PlotUtils::DefaultCVUniverse {
  public:
    // Constructor
    CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma=0);
    CVUniverse();

    // Destructor
    virtual ~CVUniverse(){};

    // All functions we write here should be 'virtual', so that the universes
    // that inherit from CVUniverse can override them.

    // ========================================================================
    // Get Weight
    //
    // We override the various weight getting functions herein in different
    // vertical systematic universe classes.
    // ========================================================================
  

  virtual double GetWeight() const;
    // double wgt_flux_and_cv=1.;
/*      double wgt_flux_and_cv=1., wgt_nrp=1., wgt_genie=1.,wgt_2p2h=1.,wgt_rpa=1.;

      // flux + cv
      std::string playlist("minervame1A");
      double Enu  = GetDouble("mc_incomingE")*1e-3;
      int nu_type = GetInt("mc_incoming");
       wgt_flux_and_cv = GetFluxAndCVWeight( Enu, nu_type);
     //  wgt_flux_and_cv = 1.0;
      // genie
      wgt_genie = GetGenieWeight();

      // mnvtune -- non-res pi
        wgt_nrp = GetNonResPiWeight();

      //2p2h
      //  wgt_2p2h = Get2p2hWeight(q0, q3);

      //RPA
        //wgt_rpa = GetRPAWeight(q0, q3);

   //   return wgt_flux_and_cv;
    return wgt_flux_and_cv*wgt_genie*wgt_nrp;

*/                                    
/*
//   const bool do_warping = true;
   double wgt_flux=1., wgt_2p2h=1.;
   double wgt_rpa=1.,   wgt_nrp=1.,  wgt_lowq2=1.;
   double wgt_genie=1., wgt_mueff=1.;
   double wgt_anisodd=1.;
 
   // genie
   wgt_genie = GetGenieWeight();
   //if (do_warping)
   //  wgt_genie = GetGenieWarpWeight();
 
   // flux
   double Enu  = GetDouble("mc_incomingE")*1e-3;
   int nu_type = GetInt("mc_incoming");
   wgt_flux = GetFluxAndCVWeight(Enu, nu_type);
 
   // 2p2h
   double q0=Getq0True(); // MeV.
   double q3=Getq3True(); // MeV.
   q0 = q0/1000.; // pass to function as GeV
   q3 = q3/1000.; // pass to function as GeV
   wgt_2p2h = Get2p2hWeight(q0, q3);
 
   // rpa
   wgt_rpa = GetRPAWeight(q0, q3);
 
   // non-res pi
   wgt_nrp = GetNonResPiWeight();
   //if (do_warping)
   //  wgt_nrp = GetNonResPiWarpWeight();
 
   // MINOS efficiency
 //  if (!m_is_truth && GetBool("isMinosMatchTrack"))
   //  wgt_mueff = GetMinosEfficiencyWeight(GetDouble("NukeCC_minos_trk_p")/1000., 
     //                                     GetThetamuDeg());
 
   //if (do_warping) {
   //  double q2 = GetQ2True();
   //  q2 = q2/1000000.; // pass to function as GeV^2
   //  wgt_lowq2 = GetLowQ2PiWarpWeight(q2, CCNuPionIncShifts::kLowQ2PiChannel); 
   //}
 
   //h_RPA_wgts->Fill(wgt_rpa);
 
   // aniso delta decay weight -- currently being used for warping
   //if (do_warping)
    // wgt_anisodd = GetVecElem("truth_genie_wgt_Theta_Delta2Npi",4);
 
   return wgt_genie*wgt_flux*wgt_2p2h*wgt_rpa*wgt_nrp*wgt_lowq2*wgt_mueff;

    }

*/
    // ========================================================================
    // Get Variable Functions
    // Write a virtual "Get" function for _any_ variable (coming directly from a
    // branch, or composed of several branches) that will be laterally shifted
    // or affected by the lateral shift of a systematic universe.
    //
    // We override some or all of these function in different systematic
    // universe classes located in LateralSystematics.h.
    // ========================================================================
    enum t_HelicityType{
      kAny = 0,
      kNeutrino,
      kAntiNeutrino
    };

//double NukeCC_leptonE[4]={1,1,1,GetVecElem("NukeCC_leptonE",3)};
//double NukeCC_leptonE_h = GetVecElem("NukeCC_leptonE",3);

    virtual double GetMuonE() const { return GetVecElem("NukeCC_leptonE",3); }

    virtual  Long64_t GetMaxEntries(); 
      virtual void     Show(Long64_t entry = -1);


virtual double calcq3(const double Q2,const double Enu,const double Emu	)const;
virtual double calcq0(const double Enu,const double Emu	)const;
      //! Turn off branches that I don't use
      //virtual bool SetBranches( bool useWeights = true );
      virtual double GetEhad() const;
          
// virtual double  GetEmu()  const ;
 virtual double GetThetamu()  const ;
 //virtual double GetQ2()     const ;
 virtual double GetQ2Reco()     const ;
 virtual double GetWReco()     const;
 virtual double GetxReco()     const;
 virtual double GetyReco()     const;

virtual double GetEnu() const;  
// for true variables
//
virtual double GetEmuTrue() const ;
virtual double GetEnuTrue() const ; 
virtual double GetEhadTrue() const;
virtual double GetThetamuTrue( )const;

 virtual double GetQ2True()     const ;
 virtual double GetWTrue()     const;
 virtual double GetyTrue()     const;
 virtual double GetxTrue()     const;
  
virtual double GetThetamuDeg()  const;
	virtual double Getq3True()const; 
	virtual double Getq0True()const; 
    //!Set up a vector of bad events
      std::vector< pair< pair<int, int>, pair<int, int> > > InitBadEventVector();

      //===== Public member variables =====//
      //! Is this tree made from MC?
      bool isMC;


virtual double calcRecoQ2(const double Enu,const double Emu,const double Thetamu)const;	
virtual double calcWReco(const double Q2,const double Ehad) const;

virtual double calcTrueQ2(const double EnuTrue,const double EmuTrue,const double ThetamuTrue)const;	
virtual double calcWTrue(const double Q2True,const double EhadTrue) const;


virtual double calcXTrue(const double Q2True,const double EnuTrue, const double EmuTrue)const;

virtual double calcXReco(const double Q2,const double Enu, const double Emu)
const ;


virtual double calcYTrue(const double EnuTrue,const double EhadTrue)const;	

virtual double calcYReco(const double Enu,const double Ehad)const;	

      //===== BEGIN ANALYSIS CUTS =====//
      virtual bool PassHelicityCut( t_HelicityType h );

       
      virtual bool PassMuCurveCut( double minCut=0 );
      
      virtual bool PassMuCoilCut( double coilR, double maxR );
      virtual bool PassMuCoilCut( );

      virtual bool PassMuQualityCut( int qual = 2 );


      //! Do the true nu and true theta mu cuts
      virtual bool PassTrueMuEnergyCut( );


   //   virtual bool PassMuEnergyCut( );
      virtual bool PassMuEnergyCut( double muonE=2. );

      virtual bool PassThetaCut( double muonTheta=17. );



     

       virtual double GetVtxEnergy(  unsigned int shift = VTX_BLOB_SHIFT );  
      virtual bool PassGoodTrackingCut();
       virtual double GetCCQERecoil( unsigned int shift = VTX_BLOB_SHIFT );
      virtual bool PassZDistCut( double lowZ = 1001., double highZ = 1001. );

      virtual bool PassDistToDivisionCut( double xySep = 25. );

      //! Is the true vertex far enough away form a division of target sections?
      virtual bool PassTrueDistToDivisionCut( double xySep  = 25. );

     virtual bool PassReco(t_HelicityType h /* = HelicityType::kNeutrino */);

      //! Check that the ccqe recoil energy is greater than this value in MeV
      virtual bool PassCCQERecoilCut( const double minE = MIN_CCQE_RECOIL_E );
      //! Check if passes the inelastic cuts using branch variables (when USE_INEL_CUT is true) 
      virtual bool PassInelasticCut();
      //! Check if passes the inelastic cuts using given variables in MeV (when USE_INEL_CUT is true) 
      virtual bool PassInelasticCut( double q2, double ccqeRecoil );

      //! Check if passes the DIS cuts using reco branch variables 
      virtual bool PassDISCut();
      
      //! Check if passes the DIS cuts using given variables in MeV 
      virtual bool PassDISCut( double q2, double W );

      

      //! Check if passes the DIS cuts using true variables in MeV 
      virtual bool PassTrueDISCut( );


      virtual bool PassTrueDISCut( double q2, double W );

      //@}
      //===== END ANALYSIS CUTS =====//


      virtual bool IsInMaterial(  int targID, int targZ, bool anyTrackerMod = true, bool isDNN = false );
      //virtual bool IsInMaterial( const MatPair& mat, bool anyTrackerMod = true );


      virtual bool IsInTrueMaterial( const int targID, const int targZ, bool anyTrackerMod = true );

      virtual double Var( const std::string& var, bool useTrue = false );

      //@}
      //===== END ANALYSIS FUNCTIONS =====//

      //! coord to U (perdendicular to section divide in tgt 1,5)
      virtual double GetU(double x, double y) const;
      //! coord transform to D (perdendicular to section divide in tgt 2)
      virtual double GetD(double x, double y) const;
      //! coord transform to C (perdendicular to carbon divide in target 3)
      virtual double GetC(double x, double y) const;

      
      virtual std::string GetXaxisTitle( const std::string& var ) const;
      virtual  std::string GetYaxisTitle( const std::string& var ) const;

       virtual double GetVarNormWidth( const std::string& var ) const;
      //! Machine Learning Stuffs
      virtual bool InHexagon( double apothem = 850. ); 
      virtual int GetTargetFromSegment( int segment, int& vtx_module, int& vtx_plane );
       
       virtual std::vector<std::string> GetStdEnergyVars( ) const;
      virtual std::vector<double> GetDISBins(const std::string& var ) const;
      virtual double GetDISVarMaxVal(const std::string& var ) const;
      virtual double GetDISVarMinVal( const std::string& var ) const;
      virtual bool IsInTargetSection( int targetID, int targetZ, double x, double y ) const;
      virtual int GetTargetPlane( int targetID ) const;
      //! Map Plane to its respective Target ID
};
}

#endif
