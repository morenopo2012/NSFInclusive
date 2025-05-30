#ifndef MNV_NUKECCUTILSNSF_cxx
#define MNV_NUKECCUTILSNSF_cxx 1

//#include "NukeCCUtilsNSF.h"
#include "../include/NukeCCUtilsNSF.h"

//#include "../include/NukeCC_Cuts.h"
//#include "CCQENuUtilsNSF.h"
#include "../include/CVUniverse.h"
//#include "include/GlobalIncludes.h" 
#include "../include/LateralSystematics.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
//#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/MinosEfficiencySystematics.h"

#include "include/CondorInput.h"

#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "MinervaUnfold/MnvResponse.h"
#include "MinervaUnfold/MnvUnfold.h"

using namespace CondorUtils;
using namespace NUKECC_ANA; 

NukeCCUtilsNSF::NukeCCUtilsNSF(string playlist){
    cvhistos1D.clear();
    cvhistos2D.clear();
     cutter = new NukeCC_Cuts();
//   GlobalParameters::Get().m_useFluxConstraint = false;
   MinervaUniverse::SetPlaylist(playlist);
  // GlobalParameters::Get().m_usePPFX1Flux = false;
   use_merged_files = false;
   
   if(neutrinoMode){
   incoming_pdg = 14;
   cout<<"Initializing NukeCCUtils in Neutrino Mode"<<endl;
   }
   else{
   incoming_pdg = -14;
   cout<<"Initializing NukeCCUtils in Antineutrino Mode"<<endl;
   
   }
   MinervaUniverse::SetAnalysisNuPDG(incoming_pdg);
  if( neutrinoMode ){
    //Use neutrino weights for general constructor
    //Vector declared in header file
    michel_weights.push_back(0.986);
    michel_weights.push_back(1.013);
  } else{
    //For antineutrino, set weights to 0. at present 
    michel_weights.push_back(1.0);
    michel_weights.push_back(1.0);
  }    

}
 


NukeCCUtilsNSF::NukeCCUtilsNSF()
{
    cvhistos1D.clear();
    cvhistos2D.clear();
     
     cutter = new NukeCC_Cuts();
 //   cutter = new CVUniverse();
  // GlobalParameters::Get().m_useFluxConstraint = false;
  // DefaultCVUniverse::SetPlaylist("minervame1A");
  // GlobalParameters::Get().m_usePPFX1Flux = false;
   use_merged_files = false;
  
    
   if(neutrinoMode){
   incoming_pdg = 14;
   cout<<"Initializing NukeCCUtils in Neutrino Mode"<<endl;
   }
   else{
   incoming_pdg = -14;
   cout<<"Initializing NukeCCUtils in Antineutrino Mode"<<endl;
   
   }
   MinervaUniverse::SetAnalysisNuPDG(incoming_pdg);
  if( neutrinoMode ){
    //Use neutrino weights for general constructor
    //Vector declared in header file
    michel_weights.push_back(0.986);
    michel_weights.push_back(1.013);
  } else{
    //For antineutrino, set weights to 0. at present 
    michel_weights.push_back(1.0);
    michel_weights.push_back(1.0);
  }    

}



 //--Destructor
 NukeCCUtilsNSF::~NukeCCUtilsNSF()
 {
    delete cutter;
     cvhistos1D.clear();
     cvhistos2D.clear();
 }



std::vector<std::string> NukeCCUtilsNSF::GetStdPlaylists( HelicityType::t_HelicityType helicity ) const
{
    std::vector<std::string> playlists;
    if( HelicityType::kAntiNeutrino != helicity )
    {
        playlists.push_back( "minerva1" );
        playlists.push_back( "minerva7" );
        playlists.push_back( "minerva9" );
        playlists.push_back( "minervame6A" );
        playlists.push_back( "minervame6B" );
        playlists.push_back( "minervame6H" );
        playlists.push_back( "minervame6I" );
        playlists.push_back( "minervame6J" );
        playlists.push_back( "minerva13C" );
                                        }
        //                                    
          if( HelicityType::kNeutrino != helicity )
          {
        playlists.push_back( "minervame1A" );

}
    
    return playlists;
}

HelicityType::t_HelicityType NukeCCUtilsNSF::GetHelicityFromPlaylist( std::string playlist ) const
{
    //if ("" == playlist)
       // playlist = fPlaylist;
    
        if(Playlist::minervame1A == playlist ||
           Playlist::minervame1B == playlist || 
           Playlist::minervame1C == playlist || 
           Playlist::minervame1D == playlist || 
           Playlist::minervame1E == playlist || 
           Playlist::minervame1F == playlist || 
           Playlist::minervame1G == playlist || 
           Playlist::minervame1M == playlist || 
           Playlist::minervame1N == playlist || 
           Playlist::minervame1O == playlist || 
           Playlist::minervame1P == playlist || 
           Playlist::minervame1L == playlist ) {
        return HelicityType::kNeutrino;
    }
       if (Playlist::minervame5A == playlist ||
        Playlist::minervame6A == playlist    ||
        Playlist::minervame6B == playlist    ||
        Playlist::minervame6H == playlist    
        ) {
        return HelicityType::kAntiNeutrino;
    }
    
    const std::vector<std::string> nus = GetStdPlaylists( HelicityType::kNeutrino );
    if ( find( nus.begin(), nus.end(), playlist ) != nus.end() )
        return HelicityType::kNeutrino;
    
    const std::vector<std::string> antinus = GetStdPlaylists( HelicityType::kAntiNeutrino );
    if ( find( antinus.begin(), antinus.end(), playlist ) != antinus.end() )
        return HelicityType::kAntiNeutrino;
    
    Warning( "NukeCCUtilsNSF::GetHelicityFromPlaylist", Form( "Playlist '%s' is not a standard nu or antinu playlist. Using any neutrino.", playlist.c_str() ) );
    return HelicityType::kAny;
}

TString NukeCCUtilsNSF::GetHistFileNamePlaylist( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ, HelicityType::t_HelicityType helicity, std::string playlist) const
{
    TString histFileName;
    
    TString variation = GetVariationTag();
    if( FileType::kAny == fType )
    {
        //histFileName += Form("/Hists_%s_t%d_z%02d_%s_%s_%s.root",
        histFileName += Form("/Hists_%s_t%d_z%02d_%s_%s.root",
                             histType.c_str(),
                             targetID, targetZ,
                             GetHelicityString( helicity ).c_str(),
                             playlist.c_str() 
                             //getenv("NUKECC_TAG"));
                             );
    }
    else
    {
        //histFileName += Form("/Hists_%s_%s_t%d_z%02d_%s_%s_%s.root",
        histFileName += Form("/Hists_%s_%s_t%d_z%02d_%s_%s.root",
                             histType.c_str(),
                             GetFileTypeString(fType).c_str(),
                             targetID, targetZ,
                             GetHelicityString( helicity ).c_str(),
                             //getenv("NUKECC_TAG"),
                             variation.Data());
                             
    }
    return histFileName; 
}

TString NukeCCUtilsNSF::GetHistFileName( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ, HelicityType::t_HelicityType helicity) const
{
    TString histFileName;
    
    TString variation = GetVariationTag();
    if( FileType::kAny == fType )
    {
        //histFileName += Form("/Hists_%s_t%d_z%02d_%s_%s_%s.root",
        histFileName += Form("/Hists_%s_t%d_z%02d_%s.root",
                             histType.c_str(),
                             targetID, targetZ,
                             GetHelicityString( helicity ).c_str() 
                             //getenv("NUKECC_TAG"));
                             );
    }
    else
    {
        //histFileName += Form("/Hists_%s_%s_t%d_z%02d_%s_%s_%s.root",
        histFileName += Form("/Hists_%s_%s_t%d_z%02d_%s_%s.root",
                             histType.c_str(),
                             GetFileTypeString(fType).c_str(),
                             targetID, targetZ,
                             GetHelicityString( helicity ).c_str(),
                             //getenv("NUKECC_TAG"),
                             variation.Data());
                             
    }
    return histFileName; 
}
/*
int UnfoldIterations( const std::string& var )
{
    int fUnfoldIter = 0; 
    if( int fUnfoldIter != - 1)
        return fUnfoldIter;
    
    if(
       "Emu" == var ||
       "Ehad" == var ||
       "Enu" == var ||
       "ThetaMu" == var
       )
        return 2;
    
    else if( "x" == var )
        return 2;
    
    return 0;
}
void SetUnfoldIter( int iter ){
    fUnfoldIter = iter;
}  
*/
TString NukeCCUtilsNSF::HistDir( bool forceDisk /* = false */, bool ignorePlaylist /*= false*/, bool needsbluearc /*= false*/ ) const
{
    TString CONDOR_DIR_HISTS = getenv("CONDOR_DIR_HISTS");
    if( ! ( CONDOR_DIR_HISTS.IsNull() || forceDisk ) )
    {
        return TString( CONDOR_DIR_HISTS + "/");
    }
    else
    {
        
        TString HISTS_ROOT;
        
        if(!needsbluearc) {
            HISTS_ROOT = getenv("HISTS");
        }
        if(needsbluearc){
            HISTS_ROOT = getenv("HISTS_BLUEARC");
        }
        cout<<" The hists directory is "<<HISTS_ROOT<<endl;
        if( HISTS_ROOT.IsNull() )
        {
            Error( "NukeCCUtilsNSF::HistDir", "Environmental variable HISTS not defined." );
            throw 1;
        }
        
        //replace minervagli with GRID_USER
        const TString minervagli = "minervagli";
        int idx = HISTS_ROOT.Index(minervagli);
        if( idx != kNPOS )
        {
            TString GRID_USER = getenv("GRID_USER");
            HISTS_ROOT.Replace(idx, idx + minervagli.Length(), GRID_USER );
            cout << "  REPLACE minervagli with GRID_USER = " << GRID_USER << ", HISTS_ROOT = " << HISTS_ROOT << endl;
        }
        
        TString NUKECC_TAG  = getenv("NUKECC_TAG");
        if( NUKECC_TAG.IsNull() )
        {
            Error( "NukeCCUtilsNSF::HistDir", "Environmental variable NUKECC_TAG not defined." );
            throw 1;
        }
        
        TString histdir(  HISTS_ROOT + "/" + NUKECC_TAG + "/"  );
       // if( fPlaylist.empty() || ignorePlaylist )
         //   histdir = TString( HISTS_ROOT + "/" + NUKECC_TAG + "/" );
        
        //cout << system( Form( "test -d %s", histdir.Data() ) ) <<endl;
        if( 0 != system( Form( "test -d %s", histdir.Data() ) ) )
        {
            int madedir = system( Form( "mkdir -m 755 -p %s", histdir.Data() ) );
            
            if( 0 != madedir )
                Error( "NukeCCUtilsNSF::HistDir", Form("Could not make hist directory '%s'", histdir.Data() ) );
        }
        
        return histdir;
    }
}


TString NukeCCUtilsNSF::GetXSecFileName() const
{
    
    //cout << fPlaylist << endl;
    if(fPlaylist == "minerva1")
        return "/minerva/data/users/joelam/NukeFiles/NukeCC_xsecMinerva1.root";
    
    else if(fPlaylist == "minerva13C")
        return "/minerva/data/users/joelam/NukeFiles/NukeCC_xsecMinerva13C.root";
    
    else{
        cout << "Playlist: " << fPlaylist << endl;
        return "/minerva/data/users/${USER}/NukeFiles/NukeCC_xsec_minervame1A_test.root";
    }
}



TString NukeCCUtilsNSF::GetVariationTag() const
{
    TString var = "";
  /* if( fSysShiftTag.empty() )
    {
        if( USE_INEL_CUT )
            var += "_inel";
        
        if( 0 != VTX_BLOB_SHIFT )
            var += Form("_VtxBlobShift%d", VTX_BLOB_SHIFT );
     */   return var;
  //  }
   // else 
  //  {
    //    var += "_" + fSysShiftTag;
  //  }
    
   // return TString( fSysShiftTag.c_str() );
}

std::string NukeCCUtilsNSF::GetHelicityString( HelicityType::t_HelicityType helicity ) const
{
    if(      HelicityType::kNeutrino     == helicity ) return "Nu";
    else if( HelicityType::kAntiNeutrino == helicity ) return "AntiNu";
    else return "AnyHelicity";
}
TString NukeCCUtilsNSF::GetHistName( const std::string& histType, FileType::t_FileType fType, const std::string& var, int targetID, int targetZ ) const
{
    TString histName;
    
    if( FileType::kAny == fType )
    {
        histName += Form( "%s_%s_t%d_z%02d",
                         histType.c_str(),
                         var.c_str(),
                         targetID,
                         targetZ
                         );
    }
    else
    {
        histName += Form( "%s_%s_%s_t%d_z%02d",
                         histType.c_str(),
                         GetFileTypeString(fType).c_str(),
                         var.c_str(),
                         targetID,
                         targetZ
                         );
    }
    return histName;
}

TString NukeCCUtilsNSF::GetHistName( const std::string& histType, FileType::t_FileType fType, const std::string& var) const
{
    TString histName;
    
    if( FileType::kAny == fType )
    {
        histName += Form( "%s_%s",
                         histType.c_str(),
                         var.c_str()
                         );
    }
    else
    {
        histName += Form( "%s_%s_%s",
                         histType.c_str(),
                         GetFileTypeString(fType).c_str(),
                         var.c_str()
                         );
    }
    return histName;
}
std::string NukeCCUtilsNSF::GetFileTypeString( FileType::t_FileType fType ) const
{
    if(      FileType::kData  == fType ) return "Data";
    else if( FileType::kMC    == fType ) return "MC";
    else if( FileType::kTruth == fType ) return "Truth";
    else if( FileType::kAny   == fType ) return "Any";
    else if( FileType::kDNNData == fType ) return "DNNData";
    else if( FileType::kDNNMC == fType ) return "DNNMC";
    else if( FileType::kDCNNData == fType ) return "DCNNData";
    else if( FileType::kDCNNMC == fType ) return "DCNNMC";
    else if( FileType::kDANNData == fType ) return "DANNData";
    else if( FileType::kDANNMC == fType ) return "DANNMC";
    else if( FileType::kNukeOnlyMC == fType ) return "NukeOnlyMC";
    else if( FileType::kNukeOnlyTruth == fType ) return "NukeOnlyTruth";
    else throw( "Unknown FileType found in GetFileTypeString" );
}

 
  void writePOT(TFile *f , double DataPOT, double MCPOT){
    //PlotUtils::ChainWrapper* chainData = util.m_data;
    //PlotUtils::ChainWrapper* chainMC = util.m_mc;
  double Data = 0;
  double MC   = 0;
  //this->getPOT(data,mc);
  TVector2 *pot = new TVector2( Data, MC );
   f->WriteTObject( pot, "pot" );
   }

///////////////For cross section stuff//////////////////////

//get average xsec of these materials
MnvH1D *NukeCCUtilsNSF::GetXSecRecoHistAvg( const std::string var, FileType::t_FileType type, MatPairs mats,  HelicityType::t_HelicityType helicity, bool isoCorrect, TFile *f1 )
{
    //get the XSec for each and average
    MnvH1D *rval(0);
    int mat = 0;
    cout<<"I'm working on getting the average cross section "<<endl;
    for( MatPairs::iterator iMat = mats.begin(); iMat != mats.end(); ++iMat )
    {
        
        cout<<"I'm on mat "<<mat<<endl;
        if( 0 == rval )
        {
            rval = GetXSecRecoHist( var, type, *iMat, helicity, isoCorrect, f1 );
            //            rval->SetBit( TH1::kIsAverage );
            cout<<"I m getting the pair "<<var<<"  "<<type<<"  "<<mat<<endl;
            rval->SetDirectory(0);
        }
        
        else
        {
            MnvH1D *tmp = GetXSecRecoHist( var, type, *iMat, helicity, isoCorrect, f1 );
            //            tmp->SetBit( TH1::kIsAverage );
            tmp->SetDirectory(0);
            /////cout<<"i'm adding the reco histogram "<<var<<"  "<<type<<"  "<<mats<<endl;
            rval->Add(tmp);
            delete tmp;
        }
        mat++;
    }
    
    //  if( rval )
    //    rval->SetBit( TH1::kIsAverage, false );
    
    ///////cout<<"I finished the cross section average making"<<var<<"  "<<type<<"  "<<mats<<endl;
    
        cout<<"Can you see me here!!! "<<mat<<endl;
    return rval;
}

MatPairs NukeCCUtilsNSF::GetStdTrackerNormPairs( int targetZ /* = 0 */ ) const
{
    MatPairs mats;
    mats.push_back( MatPair(14, 82) );
    return mats;
}

MatPairs NukeCCUtilsNSF::GetStdTrackerNormPairs( MatPair mat ) const
{
    return GetStdTrackerNormPairs( mat.first, mat.second );
}

MatPairs NukeCCUtilsNSF::GetStdTrackerNormPairs( int targetID, int targetZ ) const
{
    //get all mat pairs for this Z
    MatPairs mats = GetStdTrackerNormPairs( targetZ );
    
    //remove pairs from the wrong reference
    for( MatPairs::iterator mat = mats.begin(); mat != mats.end(); )
    {
        const int refTarg = mat->first % 10;
        if( refTarg == targetID )
            ++mat;
        else
            mat = mats.erase(mat);
    }
    
    return mats;
}


MnvH1D * NukeCCUtilsNSF::GetXSecRecoHistNorm( const std::string& var, FileType::t_FileType type, MatPair mat, HelicityType::t_HelicityType helicity, bool isoCorrect, TFile *f1 )
{
    //Get Xsec for all normalization materials and average
    MatPairs mats = GetStdTrackerNormPairs( mat );
    
    //get the XSec for each and average
    return GetXSecRecoHistAvg( var, type, mats, helicity, isoCorrect, f1 );
}

MnvH1D * NukeCCUtilsNSF::GetXSecRecoHistNormAvg( const std::string& var, FileType::t_FileType type, int i_targetZ, HelicityType::t_HelicityType helicity, bool isoCorrect, TFile *f1  )
{
    //Get Xsec for all normalization materials and average
    MatPairs mats = GetStdTrackerNormPairs( i_targetZ );
    
    //get the XSec for each and average
    return GetXSecRecoHistAvg( var, type, mats, helicity, isoCorrect, f1 );
    
}


MatPairs NukeCCUtilsNSF::GetSpecialMatPairs( int targetZ )const
{
    MatPairs mats;
    if( TRACKER_Z == targetZ )
    {
        //use all 9 faux targets (4*9=36 modules)
        for( int i = 1; i != 10; ++i )
        {
            //      if(i == 3) continue;
            int fauxID = 10*i + 4;
            mats.push_back( MatPair(fauxID,82) );
        }
    }
    else if( TRACKER_Z_ODD == targetZ )
    {
        for( int i = 1; i < 10; i+=2 )
        {
            int fauxID = 10*i + 4;
            mats.push_back( MatPair(fauxID,82) );
        }
    }
    else if( TRACKER_Z_EVEN == targetZ )
    {
        for( int i = 2; i < 10; i+=2 )
        {
            int fauxID = 10*i + 4;
            mats.push_back( MatPair(fauxID,82) );
        }
    }
    else if( TRACKER_Z_C == targetZ )
        mats = GetStdTrackerNormPairs( 6 );
    else if( TRACKER_Z_FE == targetZ )
        mats = GetStdTrackerNormPairs( 26 );
    else if( TRACKER_Z_PB == targetZ )
        mats = GetStdTrackerNormPairs( 82 );
    
    else
    {
        MatPairs tmpMats = GetStdMatPairs();
        //add in target 1
        // tmpMats.push_back( MatPair(1, 26) );
        // tmpMats.push_back( MatPair(1, 82) );
        
        for( MatPairs::iterator iMat = tmpMats.begin(); iMat != tmpMats.end(); ++iMat )
        {
            if( iMat->second == targetZ )
                mats.push_back( *iMat );
        }
    }
    cout<<"The size of my mats vector is "<<mats.size()<<endl;
    return mats;
}

MatPairs NukeCCUtilsNSF::GetStdMatPairs( int nFaux /* = 0 */ ) const
{
    MatPairs mats;
    for( int i = 0; i != nFaux + 1; ++i )
    {
        mats.push_back( MatPair( 1 + i*10, 26 ) );
        mats.push_back( MatPair( 1 + i*10, 82 ) );
        mats.push_back( MatPair( 3 + i*10,  6 ) );
        mats.push_back( MatPair( 2 + i*10, 26 ) );
        mats.push_back( MatPair( 3 + i*10, 26 ) );
        mats.push_back( MatPair( 5 + i*10, 26 ) );
        mats.push_back( MatPair( 2 + i*10, 82 ) );
        mats.push_back( MatPair( 3 + i*10, 82 ) );
        mats.push_back( MatPair( 4 + i*10, 82 ) );
        mats.push_back( MatPair( 5 + i*10, 82 ) );
    }
    return mats;
}

//sum over all targets with this Z
MnvH1D * NukeCCUtilsNSF::GetXSecRecoHistAvg( const std::string& var, FileType::t_FileType type, int targetZ, HelicityType::t_HelicityType helicity, bool isoCorrect, TFile *f1  )
{
    //if it's a set of faux targets to C, Fe, Pb use the dedicated function
    if( TRACKER_Z_C == targetZ )
        return GetXSecRecoHistNormAvg( var, type, 6, helicity, isoCorrect, f1 );
    if( TRACKER_Z_FE == targetZ )
        return GetXSecRecoHistNormAvg( var, type, 26, helicity, isoCorrect, f1 );
    if( TRACKER_Z_PB == targetZ )
        return GetXSecRecoHistNormAvg( var, type, 82, helicity, isoCorrect, f1 );
    cout<<"NukeUtils::l4219"<<endl;
    MatPairs mats = GetSpecialMatPairs(targetZ);
    //cout<<*mats.begin()<<endl;
    cout<<"Is it really true"<<endl;
    return GetXSecRecoHistAvg( var, type, mats, helicity, isoCorrect, f1 );
}
///////////////////////////////////////////////////////////


MnvH1D * NukeCCUtilsNSF::GetXSecHist( TFile *f, int targetZ, const std::string& var, HelicityType::t_HelicityType helicity, int targetID, const std::string& cutName, bool isoCorrect, bool fine )
{
    int allNucleons = 0;
    int tmpTargetZ = targetZ;
    if( 10 < targetID )
        tmpTargetZ = 0;
    TString histName = GetXSecHistName( targetID, tmpTargetZ, var, allNucleons, cutName, fine );
    cout << "getxsechist = " << histName << endl;
    std::string histNameStr = (std::string)histName;
    MnvH1D *xsec = dynamic_cast<MnvH1D*>( f->Get(histName) );
    if( xsec )
    {
        
        cout << "getevtratehist = " << histName << endl;
        
    }
    else
    {
        Warning( "NukeCCUtilsNSF::GetXSecHist", Form("Could not get xsec prediction for hist with name %s.  return null.", histName.Data() ) );
        return xsec;
    }
    
    //fold the histogram if desired
    if( targetID!=0 && !UnfoldVar(var)  )
    {
        MnvH1D *folded = new MnvH1D(*xsec);
        folded->SetName( Form( "%s_folded", xsec->GetName() ) );
        //HelicityType::t_HelicityType helicity = GetHelicityFromPlaylist(playlist);
        FoldHistoNormalized( folded, xsec, var, targetID, targetZ, helicity, false );
        delete xsec;
        xsec = folded;
    }
    
    //apply isoscalar correction after folding
    if( isoCorrect )
    {
        MnvH1D *xsecRatioIso = NukeCCUtilsNSF::GetIsoscalarCorrection( targetZ, var, targetID, cutName, fine );
        xsec->MultiplySingle( xsec, xsecRatioIso );
        xsecRatioIso->SetDirectory(0);
        delete xsecRatioIso;
    }
    
    
    return xsec;
}


TString NukeCCUtilsNSF::GetXSecHistName( int targetID, int targetZ, std::string var, int nucleon, std::string cutName, bool fine )
{
    if( "CosThetaMu" == var )
        var = "CosMu";
    
    if( fine )
        var += "_fine";
    
    std::string nucleonStr = "";
    if( PROTON_PDG == nucleon )
        nucleonStr = "_proton";
    else if( NEUTRON_PDG == nucleon )
        nucleonStr = "_neutron";
    
    //look for special tracker Z
    if( 90 <= targetZ )
        targetZ = 0;
    
    std::string nucleusStr = "";
    if( 6 == targetZ )
        nucleusStr = "carbon";
    else if( 26 == targetZ )
        nucleusStr = "iron";
    else if( 82 == targetZ )
        nucleusStr = "lead";
    else if( 0 == targetZ )
        nucleusStr = "tracker";
    else
        assert( false && "targetZ must be 0, 6, 26, 82." );
    
    if( cutName.empty() )
    {
        if( USE_INEL_CUT )
            cutName = "std_inel";
        else
            cutName = "std_dis";
    }
    
    cout<< Form("%s%s_%d_%s_%s_xsec", nucleusStr.c_str(), nucleonStr.c_str(), targetID, var.c_str(), cutName.c_str()) <<endl;
    
    return std::string( Form("%s%s_%d_%s_%s_xsec", nucleusStr.c_str(), nucleonStr.c_str(), targetID, var.c_str(), cutName.c_str()) );
    
    //    return std::string( Form("%s%s_%s_%s_xsec", nucleusStr.c_str(), nucleonStr.c_str(), var.c_str(), cutName.c_str()) );
}

///////////////////////////////////////////////////////////


MnvH1D * NukeCCUtilsNSF::GetXSecRecoHist( const std::string& var, FileType::t_FileType type, MatPair mat, HelicityType::t_HelicityType helicity, bool isoCorrect, TFile *f1  )
{
    
    TH1::AddDirectory(kFALSE);
    MnvH1D *rval(0);
    bool isMC = ( FileType::kMC == type || FileType::kDNNMC || FileType::kNukeOnlyMC );
    
    cout<<"You are not serious"<<endl; 
    
    int targetZ = mat.second;
    int targetID = mat.first;
    f1->cd();
    MnvH1D *mnv = (MnvH1D*)f1->Get( Form("unfolded_signal_%s",  var.c_str()));; 
  /* //get the unfolded and efficiency histos
    ///TString hFileName = GetHistFileName( "Unfolded", FileType::kAny, targetID, targetZ, helicity );
    //TString hFileName = Form("%s/Hists_PhysicsBackgd_without_SYS_t%d_z%02d_Nu_v1.root", outdir.c_str(), targetID, targetZ);
    //TFile *unfoldFile = TFile::Open(f1);
    
    if( ! unfoldFile )
    {
        cout << " Cannot open file: " << hFileName << endl;
        assert(false);
    }
    
    string htag = UnfoldVar(var) ? "unfolded" : "signal";
    TString hname  = GetHistName( htag, type, var, targetID, targetZ );
    cout << hname << "  " << hFileName << endl;
    
    MnvH1D* mnv = dynamic_cast<MnvH1D*>( unfoldFile->Get( hname ) );
    if( !mnv )
    {
        cout << "ERROR: could not get var " << var << " from file: " << hFileName << endl;
        assert(false);
    }
    
    hFileName = GetHistFileName( "Efficiency", FileType::kAny, targetID, targetZ, helicity );
    TFile *effFile = TFile::Open(hFileName);
    
    //efficiency correct in std generated or folded kinematics?
    string effTag = UnfoldVar(var) ? "efficiency_dis" : "efficiency_folded_dis";
    //  string effTag =  "efficiency_dis";
    hname = GetHistName( effTag, FileType::kAny, var, targetID, targetZ );
    cout<<"NukeUtils::4264: "<<hname<<endl;
    */
    //MnvH1D *eff = (MnvH1D*)effFile->Get( hname );
    MnvH1D *eff = (MnvH1D*)f1->Get( Form("unfolded_signal_%s_%s_%s",  var.c_str(), targetID, targetZ));;
    //MnvH1D *eff = (MnvH1D*)effFile->Get( hname );
    if( !eff )
    {
        cout << "ERROR: could not get efficiency for var " << var << " targetID " << targetID << ", targetZ " << targetZ << endl;
        assert(false);
    }
//    for( int bin=0;bin<mnv->GetNbinsX();bin++){
//        cout<<bin<<" plot content "<<mnv->GetBinContent(bin)<<endl;
//        cout<<bin<<" eff content  "<<eff->GetBinContent(bin)<<endl;
//    }
    
    
//    std::vector<std::string> vertNames_eff = eff->GetVertErrorBandNames();
//    std::vector<std::string> vertNames_unfolded = mnv->GetVertErrorBandNames();
//
//                for(int i=0;i<vertNames_eff.size();i++)
//                    std::cout<<"vert name eff "<<vertNames_eff[i]<<std::endl;
//                for(int i=0;i<vertNames_unfolded.size();i++)
//                    std::cout<<"vert name unfolded "<<vertNames_unfolded[i]<<std::endl;
    
    //    const double trackerAtoms= 2.22311e+27 * 92.; //C atoms used in extraction
    //  const double nucleonsInTarget = TargetUtils::Get().GetPassiveTargetNNucleons( targetID, targetZ, isMC );
    //cout << nucleonsInTarget << endl;
    //  double targetNumberScale = trackerAtoms / nucleonsInTarget; //old way. seems wrong --aen
    //     double targetNumberScale = 1.0 / nucleonsInTarget;
    
    mnv->AddMissingErrorBandsAndFillWithCV(*eff);
    eff->AddMissingErrorBandsAndFillWithCV(*mnv);
    rval =  EfficiencyCorrect( mnv, eff );
    
    std::vector<std::string> vertNames_eff_aftercorrection = eff->GetVertErrorBandNames();
    std::vector<std::string> vertNames_unfolded_aftercorrerction = mnv->GetVertErrorBandNames();
    
    //        for(int i=0;i<vertNames_eff_aftercorrection.size();i++)
    //            std::cout<<"vert name eff after "<<vertNames_eff_aftercorrection[i]<<std::endl;
    //        for(int i=0;i<vertNames_unfolded_aftercorrerction.size();i++)
    //            std::cout<<"vert name unfolded after "<<vertNames_unfolded_aftercorrerction[i]<<std::endl;
    
    //rval = mnv;
    rval->SetDirectory(0);
    rval->SetName( Form( "tmpUnfolded_%d_%d_%s_%s", targetID, targetZ, var.c_str(), isMC?"MC":"Data" ) );
    //  rval->Scale( targetNumberScale );
    cout<<"rval: "<<rval->GetBinContent(2)<<endl;
    if( isoCorrect )
    {
        MnvH1D *xsecRatioIso = GetIsoscalarCorrection( targetZ, var, targetID, "stdDIS", false );
        rval->MultiplySingle( rval, xsecRatioIso );
        xsecRatioIso->SetDirectory(0);
        delete xsecRatioIso;
    }
    
    ////delete unfoldFile;
    ////delete effFile;
    
    return rval;
}

bool NukeCCUtilsNSF::UnfoldVar( const std::string& var ) const
{
    /*if( fForceUnfold )
        return true;
    if( fForceNoUnfold )
        return false;
    */
    if( "Enu" == var )
        return true;
    
    else if( "x" == var )
        return true;
    
    else
        return false;
}

MnvH1D * NukeCCUtilsNSF::GetIsoscalarCorrection( int targetZ, std::string var, int targetID, const std::string& cutString, bool fine ) const
{
    // f = sigma_iso / sigma_sum_freenuc
    // note: these histograms have 0 stat error
    MnvH1D *isoNumerator   = GetIsoXSec( targetZ, var, targetID, cutString, fine );
    MnvH1D *isoDenominator = GetFreeNucSumXSec( targetZ, var, targetID, cutString, fine );
    
    MnvH1D *xsecRatioIso = (MnvH1D*)isoNumerator->Clone("tmpIsoRatio");
    xsecRatioIso->Divide( isoNumerator, isoDenominator );
    xsecRatioIso->SetDirectory(0); //don't delete me ROOT.  let the caller do that.
    
    isoNumerator->SetDirectory(0);
    isoDenominator->SetDirectory(0);
    delete isoNumerator;
    delete isoDenominator;
    
    return xsecRatioIso;
}


MnvH1D * NukeCCUtilsNSF::GetIsoXSec( int targetZ, std::string var, int targetID, std::string cutName, bool fine) const
{
    int tmpTargetZ = targetZ;
    if( 10 < targetID ) //Z=0 for tracker
        tmpTargetZ = 0;
    double N = 1., Z = 1.;
    GetNucleusZN(tmpTargetZ, Z, N);
    
    TFile *fRawFreeXSec = TFile::Open( GetFreeNucleonXSecFileName() );
    if(fine)
        var += "_fine";
    
    if( cutName.empty() )
    {
        
        cutName = "stdDIS";
    }
    
    //isoscalar nucleus is   A * (sig_n + sig_p)/2
    
    MnvH1D *n = (MnvH1D*)fRawFreeXSec->Get( Form("freen_%s_%s_xsec", var.c_str(), cutName.c_str()) );
    MnvH1D *p = (MnvH1D*)fRawFreeXSec->Get( Form("freep_%s_%s_xsec", var.c_str(), cutName.c_str()) );
    
    MnvH1D *rval = (MnvH1D*)p->Clone( Form("isoNucl_A%f", (N+Z) ) );
    rval->Add( n );
    //scale to per nucleon xsec
    rval->Scale( 1 / 2. );
    rval->SetDirectory(0);//don't delete me
    
    //these are normalized to bin width of 1
    //scale to desired bin width for flux integrated xsec
    if( "Enu" != var )
    {
        const double binWidthNorm = rval->GetBinWidth(1);
        rval->Scale(binWidthNorm);
    }
    //turn off any normalization to bin width
    rval->SetNormBinWidth(-1);
    
    
    n->SetDirectory(0);
    p->SetDirectory(0);
    delete n;
    delete p;
    delete fRawFreeXSec;
    
    
    //fold the histogram if desired
    if( targetID != 0 && !UnfoldVar(var) )
    {
        MnvH1D *folded = new MnvH1D(*rval);
        folded->SetName( Form( "%s_folded", rval->GetName() ) );
        ////HelicityType::t_HelicityType helicity = GetHelicityFromPlaylist(std::string playlist);
        /////FoldHistoNormalized( folded, rval, var, targetID, targetZ, helicity, false );
        // leave the gun, take the folded histogram
        delete rval;
        rval = folded;
    }
    
    //no stat errors
    for( int i = 0; i <= rval->GetNbinsX()+1; ++i )
        rval->SetBinError( i, 0. );
    
    return rval;
}


bool NukeCCUtilsNSF::GetNucleusZN( int Zin, double &Z, double &N ) const
{
    if( 6 == Zin )
    {
        Z = Z_carbon;
        N = N_carbon;
        return true;
    }
    else if( 26 == Zin )
    {
        Z = Z_iron;
        N = N_iron;
        return true;
    }
    else if( 82 == Zin )
    {
        Z = Z_lead;
        N = N_lead;
        return true;
    }
    else if( 0 == Zin || 90 < Zin )
    {
        Z = Z_scint;
        N = N_scint;
        return true;
    }
    else if( DEUTERON_Z == Zin )
    {
        Z = 1.;
        N = 1.;
        return true;
    }
    
    Error( "GetNucleusZN", Form("Cannot get Z and N for unknown Z = %d", Zin) );
    throw 1;
}



MnvH1D * NukeCCUtilsNSF::GetFreeNucSumXSec( int targetZ, std::string var, int targetID, std::string cutName, bool fine ) const
{
    int tmpTargetZ = targetZ;
    if( 10 < targetID )
        tmpTargetZ = 0;
    double N = 1., Z = 1.;
    GetNucleusZN(tmpTargetZ, Z, N);
    
    TFile *fRawFreeXSec = TFile::Open( GetFreeNucleonXSecFileName() );
    if(fine)
        var += "_fine";
    
    if( cutName.empty() )
    {
        cutName = "stdDIS";
    }
    
    
    MnvH1D *n = (MnvH1D*)fRawFreeXSec->Get( Form("freen_%s_%s_xsec", var.c_str(), cutName.c_str()) );
    MnvH1D *p = (MnvH1D*)fRawFreeXSec->Get( Form("freep_%s_%s_xsec", var.c_str(), cutName.c_str()) );
    
    
    MnvH1D *rval = (MnvH1D*)p->Clone( Form("freeNuclSum_Z%f_N%f", Z, N ) );
    rval->Reset();
    rval->Add( n, N );
    rval->Add( p, Z );
    rval->SetDirectory(0);//don't delete me
    
    //these are normalized to bin width of 1
    //scale to desired bin width for flux integrated xsec
    if( "Enu" != var )
    {
        const double binWidthNorm = rval->GetBinWidth(1);
        rval->Scale(binWidthNorm);
    }
    //turn off any normalization to bin width
    rval->SetNormBinWidth(-1);
    
    //scale to per nucleon
    rval->Scale( 1. / (N+Z) );
    
    
    n->SetDirectory(0);
    p->SetDirectory(0);
    delete n;
    delete p;
    delete fRawFreeXSec;
    
    //fold if appropriate
    if( targetID!=0 && !UnfoldVar(var)  )
    {
        MnvH1D *folded = new MnvH1D(*rval);
        folded->SetName( Form( "%s_folded", rval->GetName() ) );
        ////HelicityType::t_HelicityType helicity = GetHelicityFromPlaylist(playlist);
        ////FoldHistoNormalized( folded, rval, var, targetID, targetZ, helicity, false );
        // leave the gun, take the folded histogram
        delete rval;
        rval = folded;
    }
    
    //no stat errors
    for( int i = 0; i <= rval->GetNbinsX()+1; ++i )
        rval->SetBinError( i, 0. );
    
    return rval;
}

bool NukeCCUtilsNSF::FoldHisto(  MnvH1D* hFolded,
                          const MnvH1D *hGenerated,
                          const std::string& var, int targetID, int targetZ, HelicityType::t_HelicityType helicity, bool addSys, bool isDIS /*= true*/ ) const
{
    hFolded->Reset();
    hFolded->ClearAllErrorBands();
    
    TString histFileName = GetHistFileName( "Migration", FileType::kAny, targetID, targetZ, helicity );
    TFile fMigration( histFileName, "READ" );
    
    TString hNameResp;
    
    if(isDIS)
        hNameResp = GetHistName( "migration", FileType::kAny, var, targetID, targetZ );
    else
        hNameResp = GetHistName( "migration_inc", FileType::kAny, var, targetID, targetZ );
    
    //cout << "NukeUtils::FoldHisto Histo name for folding: " << hNameResp << endl;
    
    MnvH2D *hResp = dynamic_cast<MnvH2D*>( fMigration.Get( hNameResp ) );
    if( 0 == hResp )
    {
        Error( "FoldHisto", Form("Could not find response migration histogram for var %s, targetID %d, targetZ %d", var.c_str(), targetID, targetZ ) );
        return false;
    }
    
    //add BG scale error to migration if it exists on folded hist
   /////// if( hFolded->HasVertErrorBand( ErrorName::PlasticBG ) && ! hResp->HasVertErrorBand( ErrorName::PlasticBG ) )
   //////     hResp->AddVertErrorBandAndFillWithCV( ErrorName::PlasticBG, 2 );
    /*std::vector<string> latErrors = hResp->GetLatErrorBandNames();
     std::vector<string> vertErrors = hResp->GetVertErrorBandNames();
     for( uint i = 0; i < latErrors.size(); i++ ){
     cout << "name of errors migration  = " << latErrors[i] << endl;
     }
     for( uint i = 0; i < vertErrors.size(); i++ ){
     cout << "name of vert errors migration  = " << vertErrors[i] << endl;
     }*/
    return MinervaUnfold::MnvUnfold::Get().FoldHisto( hFolded, hResp, hGenerated, addSys );
}
bool NukeCCUtilsNSF::FoldHistoNormalized(MnvH1D* hFolded,
                                    const MnvH1D *hGenerated,
                                    const std::string& var, int targetID, int targetZ, HelicityType::t_HelicityType helicity,  bool addSys ) const
{
    MnvH1D *tmp = (MnvH1D*)hGenerated->Clone( Form( "%s_tmpFold",hGenerated->GetName()) );
    tmp->SetDirectory(0);
    const double binWidthNorm = tmp->GetBinWidth(1);
    //undo bin width normalization or else folding doesn't work
    for( int iBin = 0; iBin <= tmp->GetNbinsX()+1; ++iBin )
    {
        const double binWidth = tmp->GetBinWidth(iBin);
        const double content  = tmp->GetBinContent(iBin);
        const double err      = tmp->GetBinError(iBin);
        tmp->SetBinContent( iBin, content * binWidth / binWidthNorm );
        tmp->SetBinError( iBin,   err * binWidth / binWidthNorm );
    }
    
    bool foldOK = FoldHisto( hFolded, tmp, var, targetID, targetZ, helicity, addSys );
    
    //redo normalization
    hFolded->Scale( binWidthNorm, "width" );
    delete tmp;
    
    return foldOK;
}

TString NukeCCUtilsNSF::GetFreeNucleonXSecFileName() const
{
    const std::string FILES( getenv("FILES") );
    
    //////return CondorInput::Get().FetchInput(Form("%s/FreeNucleonIso_xsec.root",FILES.c_str()));
}

MnvH1D * NukeCCUtilsNSF::EfficiencyCorrect( MnvH1D *h, MnvH1D *eff ) const
{
    
    //add dummy universes for bg scale error
    //if( !h->HasVertErrorBand( ErrorName::BG_corr ) )
  //  h->AddVertErrorBandAndFillWithCV( ErrorName::BG_corr, 1 );
  //if( !eff->HasVertErrorBand( ErrorName::BG_corr ) )
  //    eff->AddVertErrorBandAndFillWithCV( ErrorName::BG_corr, 1 );
    
    //add dummy hist for bg stat error
    //if( !h->HasUncorrError( ErrorName::BG_uncorr ) )
    //  h->AddUncorrErrorAndFillWithCV( ErrorName::BG_uncorr );
    // if( !eff->HasUncorrError( ErrorName::BG_uncorr ) )
    //  eff->AddUncorrErrorAndFillWithCV( ErrorName::BG_uncorr );
    
    //add dummy hist of Efficiency stat error
    //if( !h->HasUncorrError( ErrorName::EffStat ) ){
    //  cout<<"there isn't an eff stat error on the histogram "<<endl;
    //  h->AddUncorrErrorAndFillWithCV( ErrorName::EffStat );
        
    //  cout<<"Does the efficiency have an eff stat error? "<<eff->HasUncorrError( ErrorName::EffStat ) <<endl;
    //    }
    
    
    
    //add dummy hist of Purity stat error
    //if( !eff->HasUncorrError( ErrorName::PurStat ) )
  //  eff->AddUncorrErrorAndFillWithCV( ErrorName::PurStat );
    eff->AddMissingErrorBandsAndFillWithCV(*h);
  
    MnvH1D *effCorrected = dynamic_cast<MnvH1D*>( h->Clone( Form( "%s", h->GetName() ) ) );
    effCorrected->Divide( h, eff );
    
    return effCorrected; //failures?
}

bool NukeCCUtilsNSF::DivideByFlux( MnvH1D* h, std::string playlist, bool reweightFlux ) const
{
    
    //get helicity from playlist
    //playlist = GetStdPlaylists();
    HelicityType::t_HelicityType helicity = GetHelicityFromPlaylist(playlist);
    
    //get incoming pdg from helicity
    int incoming_pdg = 14;
    if( helicity == HelicityType::kAntiNeutrino ) incoming_pdg = -14;
    
    //Initialize flux reweighter here
    FluxReweighter* fluxReweighter = new FluxReweighter( 14, USE_NUE, playlist, FluxReweighter::gen2thin, FluxReweighter::g4numiv6, N_VERT_UNIVERSES );
    
    MnvH1D* hFlux;
    
    if( reweightFlux ) hFlux = fluxReweighter->GetFluxReweighted(incoming_pdg);
    else hFlux = fluxReweighter->GetFluxGenerated(incoming_pdg);
    double mcpot = 1.0;
    
    MnvH1D* hFluxrebin = fluxReweighter->GetRebinnedFluxReweighted(incoming_pdg, h);//, reweightFlux );
    
    //hFluxrebin->SetNormBinWidth(-1.);
    if( reweightFlux ) hFluxrebin->AddMissingErrorBandsAndFillWithCV(*h);
    double scaleFactor = 1E-4; //need to convert from m^2 to cm^2
    //    double target = 2.2231E27*92;
    //    double mcPOT = GetMCPOT(playlist);
    
    hFlux->Scale(scaleFactor);
    hFluxrebin->Scale(scaleFactor);
    
    hFluxrebin->SaveAs("Rebinned_Flux.root");
    
    //hFluxDraw = dynamic_cast<MnvH1D*>( hFluxrebin->Clone("hFluxDraw") );
    //hFluxDraw->Scale(target);
    //hFluxDraw->Scale(mcPOT);
    
    //    cout<<"Hey, I scaled my flux! "<<endl;
    
    cout<<"My flux histogram has "<<hFluxrebin->GetNbinsX()<<" bins "<<endl;
    cout<<"My numerator histogram has "<<h->GetNbinsX()<<" bins "<<endl;
    h->Divide(h,hFluxrebin);
    
    return true;
}

bool NukeCCUtilsNSF::DivideByIntegratedFlux( MnvH1D* h, std::string playlist, bool reweightFlux ) const
{
    
    //get helicity from playlist
    //NukeUtils::Get().fPlaylist = playlist;
    HelicityType::t_HelicityType helicity = GetHelicityFromPlaylist(playlist);
    
    //get incoming pdg from helicity
    int incoming_pdg = 14;
    if( helicity == HelicityType::kAntiNeutrino ) incoming_pdg = -14;
    
    
    FluxReweighter* fluxReweighter = new FluxReweighter( 14, USE_NUE, playlist, FluxReweighter::gen2thin, FluxReweighter::g4numiv6, N_VERT_UNIVERSES );
    
    MnvH1D* hFlux, *hFluxrebin;
    
    const double minE = USE_INEL_CUT ?  MIN_NUE_INEL_CUT : MIN_RECO_E;
    //hFluxrebin = frw->GetIntegratedFluxReweighted( incoming_pdg, h, minE/1000., MAX_RECO_E/1000., reweightFlux );
    hFluxrebin = fluxReweighter->GetIntegratedFluxReweighted( incoming_pdg, h, 0, 200.);//, reweightFlux );
    hFluxrebin->SaveAs("Integrated_Flux.root");
    if( reweightFlux ) hFluxrebin->AddMissingErrorBandsAndFillWithCV(*h);
    double scaleFactor = 1E-4; //need to convert from m^2 to cm^2
    
    
    hFluxrebin->Scale(scaleFactor);
    
    h->Divide(h,hFluxrebin);
    
    return true;
}

double NukeCCUtilsNSF::GetTotalScatteringCenters( int targetZ, bool isMC ) const
{
    double Nucleons = -999.;
    
    if(targetZ==6)
        Nucleons = TargetUtils::Get().GetPassiveTargetNNucleons( 3, targetZ, isMC );
    if(targetZ==26)
        Nucleons = TargetUtils::Get().GetPassiveTargetNNucleons( 1, targetZ, isMC )
        + TargetUtils::Get().GetPassiveTargetNNucleons( 2, targetZ, isMC )
        + TargetUtils::Get().GetPassiveTargetNNucleons( 3, targetZ, isMC )
        + TargetUtils::Get().GetPassiveTargetNNucleons( 5, targetZ, isMC );
    if(targetZ==82)
        Nucleons =  TargetUtils::Get().GetPassiveTargetNNucleons( 1, targetZ, isMC )
        + TargetUtils::Get().GetPassiveTargetNNucleons( 2, targetZ, isMC )
        + TargetUtils::Get().GetPassiveTargetNNucleons( 3, targetZ, isMC )
        + TargetUtils::Get().GetPassiveTargetNNucleons( 4, targetZ, isMC )
        + TargetUtils::Get().GetPassiveTargetNNucleons( 5, targetZ, isMC );
    if(targetZ>90 )
        Nucleons = TargetUtils::Get().GetTrackerNNucleons( 12.0, isMC);
        //Nucleons = TargetUtils::Get().GetTrackerNNucleons( 108.0, isMC);
    //        Nucleons = TargetUtils::Get().GetTrackerNNucleons( 77.0, isMC);
    
    return Nucleons;
    
    
}

std::string NukeCCUtilsNSF::GetTotalXSecString(int targetID, int targetZ ) const
{
    string z = "A";
    if( targetID > 10 )       z = "CH";
    else if( 90 <= targetZ )  z = "CH";
    else if( targetZ == 6  )  z = "C";
    else if( targetZ == 26 )  z = "Fe";
    else if( targetZ == 82 )  z = "Pb";
    
    return string( Form("#sigma^{%s}", z.c_str() ) );
}

std::string NukeCCUtilsNSF::GetTotalXSecUnits( ) const
{
    return string(" ( 10^{-38} cm^{2} / Nucleon )");
}
 
std::string NukeCCUtilsNSF::GetXSecString( const std::string& var, const MatPair& mat ) const
{
    return GetXSecString( var, mat.first, mat.second );
}

std::string NukeCCUtilsNSF::GetXSecString( const std::string& var, const int targetID, const int targetZ ) const
{
    string z = "A";
    if( targetID > 10 )       z = "CH";
    else if( 90 <= targetZ )  z = "CH";
    else if( targetZ == 6  )  z = "C";
    else if( targetZ == 26 )  z = "Fe";
    else if( targetZ == 82 )  z = "Pb";
    
    string super = z;
    
    string denom = "";
    if( var == "Emu" )           denom = "E_{#mu}";
    else if( var == "Enu" )      denom = "E_{#nu}";
    else if( var == "Ehad" )     denom = "E_{had}";
    else if( var == "CCQE-Recoil" )     denom = "E_{recoil}^{ccqe}";
    else if( var == "VtxEnergy" )       denom = "E_{vertex}";
    else if( var == "Emum" )     denom = "E_{mum}";
    else if( var == "ETheta" )   denom = "E_{#mu}*(1-cos(#theta_{#mu}))";
    else if( var == "Q2" )       denom = "Q^{2}";
    else if( var == "W" )        denom = "W";
    else if( var == "y" )        denom = "y";
    else if( var == "x" )        denom = "x";
    else if( var == "PhiMu" )    denom = "#phi_{#mu}";
    else if( var == "CosThetaMu" )  denom = "cos(#theta_{#mu})";
    else if( var == "ThetaMu" )  denom = "#theta_{#mu}";
    else if( var == "ThetaMuX" ) denom = "#theta_{#mu,x}";
    else if( var == "ThetaMuY" ) denom = "#theta_{#mu,y}";
    else denom = "X";
    
    const string numer = Form("d#sigma^{%s}", super.c_str() );
    
    string rval = Form( "#frac{%s}{d%s}", numer.c_str(), denom.c_str() );
    return rval;
    
}


//==============================================================
// book histograms
//==============================================================


#endif
