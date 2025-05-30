#ifndef MNV_NUKECCUTILSNSF_cxx
#define MNV_NUKECCUTILSNSF_cxx 1

//#include "NukeCCUtilsNSF.h"
#include "include/NukeCCUtilsNSF.h"

//#include "../include/NukeCC_Cuts.h"
//#include "CCQENuUtilsNSF.h"
#include "include/CVUniverse.h"
#include "include/GlobalIncludes.h" 
#include "include/LateralSystematics.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"

#include "PlotUtils/MinosEfficiencySystematics.h"


#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "MinervaUnfold/MnvResponse.h"
#include "MinervaUnfold/MnvUnfold.h"

using namespace NUKECC_ANA; 

//--Default constructor
/*
CCQENuUtilsNSF::CCQENuUtilsNSF(  bool useFluxConstraint,string playlist )  
{ 
   cvhistos1D.clear();
   cvhistos2D.clear();
   
   //cutter = new CVUniverse();
  // GlobalParameters::Get().m_useFluxConstraint = useFluxConstraint;
  DefaultCVUniverse::SetNuEConstraint(useFluxConstraint);
  DefaultCVUniverse::SetPlaylist(playlist);
  // GlobalParameters::Get().m_usePPFX1Flux = false;
  // use_merged_files = use_mergedFiles;
   
   if(neutrinoMode){

   incoming_pdg = 14;
   cout<<"Initializing CCQENuUtils in Neutrino Mode"<<endl;
   }
   else{
   incoming_pdg = -14;
   cout<<"Initializing CCQENuUtils in Antineutrino Mode"<<endl;
   
   }
  DefaultCVUniverse::SetAnalysisNuPDG(incoming_pdg);
  
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
*/
 
NukeCCUtilsNSF::NukeCCUtilsNSF(string playlist){
    cvhistos1D.clear();
    cvhistos2D.clear();
     cutter = new NukeCC_Cuts();
//   GlobalParameters::Get().m_useFluxConstraint = false;
   DefaultCVUniverse::SetPlaylist(playlist);
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
   DefaultCVUniverse::SetAnalysisNuPDG(incoming_pdg);
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
   DefaultCVUniverse::SetAnalysisNuPDG(incoming_pdg);
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
 //   if ("" == playlist)
 //       playlist = fPlaylist;
    
        if(Playlist::minervame1A == playlist ||
           Playlist::minervame1B == playlist ||
           Playlist::minervame1L == playlist ) {
        return HelicityType::kNeutrino;
    }
       if (Playlist::minervame5A == playlist ||
        Playlist::minervame6A == playlist
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

TString NukeCCUtilsNSF::GetHistFileName( const std::string& histType, FileType::t_FileType fType, int targetID, int targetZ, HelicityType::t_HelicityType helicity) const
{
    TString histFileName;
    
    TString variation = GetVariationTag();
    if( FileType::kAny == fType )
    {
        histFileName += Form("/Hists_%s_t%d_z%02d_%s_%s_%s.root",
                             histType.c_str(),
                             targetID, targetZ,
                             GetHelicityString( helicity ).c_str(),
                             getenv("NUKECC_TAG"),variation.Data());
    }
    else
    {
        histFileName += Form("/Hists_%s_%s_t%d_z%02d_%s_%s_%s.root",
                             histType.c_str(),
                             GetFileTypeString(fType).c_str(),
                             targetID, targetZ,
                             GetHelicityString( helicity ).c_str(),
                             getenv("NUKECC_TAG"),variation.Data());
                             
    }
    return histFileName; 
}


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

 //============================================================================
 // bookHistosCV
 //===========================================================================
 void NukeCCUtilsNSF::bookHistosCV( HistWrapper<CVUniverse>** h, string var_name, string title, std::vector<double> nbin, double min, const double binWidth, std::map< std::string, std::vector<CVUniverse*> > error_bands ){
   MnvH1D *hist[nHistosMy];
   //if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
     for( unsigned int i = 0; i < nHistosMy; ++i ){
       hist[i] = new MnvH1D( ( var_name + '_' + namesMy[i] ). c_str(),title.c_str(),
			  nbin.size()-1, min ,binWidth); 
       h[i] = new HistWrapper<CVUniverse>( hist[i],error_bands );
       
     }
  // }
  // else{
    // for( unsigned int i = 0; i < nHistosInc; ++i )
      // h[i] = new HistWrapper<CVUniverse>( hist[i],error_bands );
  // }
    cvhistos1D[var_name] = h;
 }
 void NukeCCUtilsNSF::bookHistosCV( HistWrapper<CVUniverse>** h, string var_name, string title, std::vector<double> nbin, double min, double max, const double binWidth, std::map< std::string, std::vector<CVUniverse*> > error_bands ){
   MnvH1D *hist[nHistosMy];
   //if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
     for( unsigned int i = 0; i < nHistosMy; ++i ){
       hist[i] = new MnvH1D( ( var_name + '_' + namesMy[i] ). c_str(),title.c_str(),
			  nbin.size()-1, min , max, binWidth); 
       h[i] = new HistWrapper<CVUniverse>( hist[i],error_bands );
       
     }
  // }
  // else{
    // for( unsigned int i = 0; i < nHistosInc; ++i )
      // h[i] = new HistWrapper<CVUniverse>( hist[i],error_bands );
  // }
    cvhistos1D[var_name] = h;
 }
void NukeCCUtilsNSF::bookHistosCV( HistWrapper<CVUniverse>** h, string var_name, string title, std::vector<double> nbin, std::map< std::string, std::vector<CVUniverse*> > error_bands ){
   MnvH1D *hist[nHistosMy];
   //if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
   for( unsigned int i = 0; i < nHistosMy; ++i ){
      hist[i] = new MnvH1D( ( var_name + '_' + namesMy[i] ). c_str(),title.c_str(),
                                                  nbin.size()-1,&nbin[0]);
      h[i] = new HistWrapper<CVUniverse>( hist[i],error_bands );
     }
   // }
   // else{
   // for( unsigned int i = 0; i < nHistosInc; ++i )
   // h[i] = new HistWrapper<CVUniverse>( hist[i],error_bands );
   // }
   cvhistos1D[var_name] = h;
  }


/////////////////////////////////Amit style//////////////////////
 void NukeCCUtilsNSF::bookHistosCV( HistWrapper<CVUniverse>** h, string var_name, string title, axis_binning xbins, std::map< std::string, std::vector<CVUniverse*> > error_bands ){
  MnvH1D *hist[nHistosMy];
  // if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
     for( unsigned int i = 0; i < nHistosMy; ++i ){
       
       hist[i] = new  MnvH1D( ( var_name + '_' + names[i] ).c_str(), title.c_str(),
			   xbins.nbins, &( xbins.bin_edges[0] ), xbins.default_width );
			   
       h[i] = new HistWrapper<CVUniverse>(hist[i],error_bands);
	   
	}		   
  // }
   //else{
   //  for( unsigned int i = 0; i < nHistosInc; ++i ){
     //  hist[i] = new  MnvH1D( ( var_name + '_' + namesInc[i] ).c_str(), title.c_str(),
	//		   xbins.nbins, &( xbins.bin_edges[0] ), xbins.default_width);
       //h[i] = new HistWrapper<CVUniverse>(hist[i],error_bands)	;   
			   
    //  }
   //}
    cvhistos1D[var_name] = h;
 }
////////////////////////////////////////////////////////////////
  void NukeCCUtilsNSF::bookHistosCV( Hist2DWrapper<CVUniverse>** h, string var_name, string x_title, std::vector<double> xnbin, double xmin, double xmax, std::vector<double> ynbin, double ymin, double ymax, std::map< std::string, std::vector<CVUniverse*> > error_bands ){
   MnvH2D *hist[nHistosMy]; 
  // if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
     for( unsigned int i = 0; i < nHistosMy; ++i ){
       hist[i] = new MnvH2D( ( var_name + '_' + namesMy[i] ).c_str(),x_title.c_str(),
                            xnbin.size()-1,xmin, xmax, ynbin.size()-1,ymin, ymax);
      h[i] = new Hist2DWrapper<CVUniverse>(hist[i],error_bands);	
      
      		 
     }
    cvhistos2D[var_name] = h;
 }
  void NukeCCUtilsNSF::bookHistosCV( Hist2DWrapper<CVUniverse>** h, string var_name, string x_title, std::vector<double> xnbin, std::vector<double> ynbin, std::map< std::string, std::vector<CVUniverse*> > error_bands ){
   MnvH2D *hist[nHistosMy]; 
  // if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
     for( unsigned int i = 0; i < nHistosMy; ++i ){
       hist[i] = new MnvH2D( ( var_name + '_' + namesMy[i] ).c_str(),x_title.c_str(),
                            xnbin.size()-1, &xnbin[0], ynbin.size()-1, &ynbin[0]);
      h[i] = new Hist2DWrapper<CVUniverse>(hist[i],error_bands);	
      
      		 
     }
    cvhistos2D[var_name] = h;
 }
///////////////////////////////////////Amit style//////////////////////////// 
 void NukeCCUtilsNSF::bookHistosCV( Hist2DWrapper<CVUniverse>** h,string var_name,string title,axis_binning xbins, axis_binning ybins, std::map< std::string, std::vector<CVUniverse*> > error_bands ){
   Double_t xbinsarray[xbins.nbins];
   MnvH2D *hist[nHistosMy];
   for( unsigned int i = 0; i <= xbins.nbins; i++ )
     xbinsarray[i] = xbins.bin_edges[i];
   Double_t ybinsarray[ybins.nbins];
   for( unsigned int i = 0; i <= ybins.nbins; i++ )
     ybinsarray[i] = ybins.bin_edges[i];
//   if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
     for( unsigned int i = 0; i < nHistosMy; i++ ){
      hist[i] = new MnvH2D((var_name+Form("_%s",names[i].c_str())).c_str(),title.c_str(),xbins.nbins,xbinsarray,ybins.nbins,ybinsarray);
       
       h[i] = new Hist2DWrapper<CVUniverse>(hist[i],error_bands); 
       
     }  
  // }
 
  //else{
    // for( unsigned int i = 0; i < nHistosInc; i++ ){
     //  hist[i] = new MnvH2D((var_name+Form("_%s",namesInc[i].c_str())).c_str(),title.c_str(),xbins.nbins,xbinsarray,ybins.nbins,ybinsarray);
     //  h[i] = new Hist2DWrapper<CVUniverse>(hist[i],error_bands);
 //  }
     //}
    cvhistos2D[var_name] = h;
 }








/////////////////////////////////////////////////////////////////////////
//I will only have the fill for 1D and 2D histograms for now....

 void NukeCCUtilsNSF::fillHistos( HistWrapper<CVUniverse>**h, double var, CVUniverse* universe, double w ){  // fill mc
  // if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
 //   NukeCC_Cuts* cutter = new NukeCC_Cuts();
//cout<<"Huma I ams here"<<endl;

        h[kIncMCMy]->univHist(universe)->Fill( var, w );

       if(universe->ShortName()=="cv"){ //only fill the CV histogram
       if( cutter->PassDISCut( universe ) )  h[kIncDISMy]->univHist(universe)->Fill( var, w );
       if( cutter->PassLowWCut( universe ) )  h[kInclowWMy]->univHist(universe)->Fill( var, w );
       if( cutter->PassLowQ2Cut( universe ) ) h[kInclowQ2My]->univHist(universe)->Fill( var, w );
       if( cutter->PassLowQ2Trans( universe ) ) h[kInclowQ2transMy]->univHist(universe)->Fill( var, w );
       if( cutter->Passtrans( universe ) ) h[kInctransMy]->univHist(universe)->Fill( var, w );
//cout<<"Huma I ams here 4"<<endl;   
   //  if( !cutter->passTrueCCQE( universe ) && !cutter->passTrueCCRES( universe ) && !cutter->passTrueCCDIS( universe ) && !cutter->passTrueMEC( universe ) )
     //  h[kOTH]->univHist(universe)->Fill( var, w);
     
//cout<<"Huma I ams here 6"<<endl;   
   }

 }


void NukeCCUtilsNSF::fillHistos( Hist2DWrapper<CVUniverse>**h, double var_x, double var_y,CVUniverse* universe, double w ){  // fill mc
  //if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
    //std::cout<<"Filling histograms "<<var_x<<" "<<var_y<<" "<<w<<std::endl;
  //  h[kMC]->univHist(universe)->Fill( var_x, var_y, w );

        h[kIncMCMy]->univHist(universe)->Fill( var_x, var_y, w );
      if(universe->ShortName()=="cv"){ //only fill the CV histogram


       if( cutter->PassDISCut( universe ) )  h[kIncDISMy]->univHist(universe)->Fill( var_x, var_y, w );

     if( cutter->PassLowWCut( universe ) )  h[kInclowWMy]->univHist(universe)->Fill( var_x, var_y, w );

//cout<<"Huma I ams here 4"<<endl;   
   //  if( !cutter->passTrueCCQE( universe ) && !cutter->passTrueCCRES( universe ) && !cutter->passTrueCCDIS( universe ) && !cutter->passTrueMEC( universe ) )
     //  h[kOTH]->univHist(universe)->Fill( var, w);
     
      if( cutter->PassLowQ2Cut( universe ) ) h[kInclowQ2My]->univHist(universe)->Fill( var_x, var_y, w );
     
      if( cutter->PassLowQ2Trans( universe ) ) h[kInclowQ2transMy]->univHist(universe)->Fill( var_x, var_y, w );

      if( cutter->Passtrans( universe ) ) h[kInctransMy]->univHist(universe)->Fill( var_x, var_y, w );

 } 
}

void NukeCCUtilsNSF::syncCVHistos(HistWrapper<CVUniverse>**h){
    //the sync hist function....
    //for(unsigned int i=kIncMCMy;i<nHistosMy;++i)h[i]->SyncCVHistos();
    for(unsigned int i=kData;i<nHistosMy;++i)h[i]->SyncCVHistos();
}

//void NukeCCUtilsNSF::syncCVHistos(HistWrapper<CVUniverse>*h){
//    //the sync hist function....
//    h->SyncCVHistos();
//}

void NukeCCUtilsNSF::syncCVHistos(Hist2DWrapper<CVUniverse>**h){
    //the sync hist function....
    for(unsigned int i=kData;i<nHistosMy;++i)h[i]->SyncCVHistos();
    //for(unsigned int i=kIncMCMy;i<nHistosMy;++i)h[i]->SyncCVHistos();
}

//void NukeCCUtilsNSF::syncCVHistos(Hist2DWrapper<CVUniverse>*h){
    //the sync hist function....
//    h->SyncCVHistos();
//}
//
//
//
/*
//setup a 2D MnvRepsonse object....
//template<class T>
void NukeCCUtilsNSF::setupResponse( MinervaUnfold::MnvResponse*& response, string var_name, string title,std::vector<double> x_reco, std::vector<double> y_reco,std::vector<double> x_truth, std::vector<double> y_truth, std::map<const std::string,const int>systematics ){
   response = new MinervaUnfold::MnvResponse((var_name+"_qelike" ).c_str(), title.c_str(), x_reco, y_reco, x_truth, y_truth,systematics);
   cout << response << endl;

}

//fill the response object...
//template<class T>

void NukeCCUtilsNSF::fillResponse( MinervaUnfold::MnvResponse* response,CVUniverse *cv, double var_x, double var_y, double var_xtrue, double var_ytrue,  double w, const std::string error_name,const int univ  ){  // fill mc

  if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
    if( cutter->passTrueCCQELike( cv ) ) response->Fill( var_x, var_y, var_xtrue, var_ytrue,error_name, univ, w );
  }
  else{
    if(cutter->passTrueCC( cv ))  response->Fill( var_x, var_y, var_xtrue, var_ytrue, error_name, univ, w );
  }
}


//functions for the MigrationMatrices......

 //==============================================================
 // Setup Response
 //==============================================================



void NukeCCUtilsNSF::getResponseObjects(MinervaUnfold::MnvResponse* response, MnvH2D*& h_migration, MnvH2D*& h_reco, MnvH2D*& h_truth)
{
  bool status = false;
  status = response->GetMigrationObjects( h_migration, h_reco, h_truth );
  if (!status){
    if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive) Error("CCQENuUtils::getResponseObjects", Form("Could not get migration objects for %s histograms", names[kQELike].c_str() ) ); 
    else  Error("CCQENuUtils::getResponseObjects", Form("Could not get migration objects for %s histograms", names[kIncCC].c_str() ) ); 
  }



}
 
*/

void NukeCCUtilsNSF::fillHistos( HistWrapper<CVUniverse>** h, double var ){  // fill data
 // if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
    h[kData]->hist->Fill( var );
  }

void NukeCCUtilsNSF::fillHistos( Hist2DWrapper<CVUniverse>** h, double var_x, double var_y ){  // fill data
  //if(GlobalParameters::Get().m_analysisType!= kInclusive && GlobalParameters::Get().m_analysisType != k3DInclusive){
    h[kData]->hist->Fill( var_x,var_y );
  }
 
ChainWrapper* NukeCCUtilsNSF::GetChainWrapperMCPointer(string playlist,string tree_name){
  std::string mybase = getenv("MY_NSFNUKECC");
  std::string mc_playlist_path = mybase+"/include/playlists/NukeCC_"+playlist+"_MC_Inextinguishable_merged.txt";
  ChainWrapper* chain = makeChainWrapperPtr(mc_playlist_path,tree_name); 
 return chain;
} 
 
 
ChainWrapper* NukeCCUtilsNSF::GetChainWrapperDataPointer(string playlist,string tree_name){
  std::string mybase = getenv("MY_NSFNUKECC");
  std::string data_playlist_path = mybase+"/include/playlists/NukeCC_"+playlist+"_DATA_Inextinguishable_merged.txt";
  ChainWrapper* chain = makeChainWrapperPtr(data_playlist_path, tree_name);
  return chain;
} 


void NukeCCUtilsNSF::getPOT(double& total_pot_data,double& total_pot_mc){
     ChainWrapper *chain_mc =  GetChainWrapperMCPointer(DefaultCVUniverse::GetPlaylist(),"Meta");
     ChainWrapper *chain_data = GetChainWrapperDataPointer(DefaultCVUniverse::GetPlaylist(),"Meta");
     total_pot_data = setPOTData(chain_data);
     total_pot_mc = setPOTMC(chain_mc);
}

double NukeCCUtilsNSF::setPOTData(ChainWrapper* _chain){
   TChain *chain = _chain->GetChain();
   double tot_pot_data = 0.0;
   double pot_data = 0.0;
   TBranch *b_pot_data;
   chain->SetBranchAddress("POT_Used",&pot_data,&b_pot_data);
   Int_t tot_entries = chain->GetEntries();
   double global_pot_data = 0.0;
   double global_pot=0.0;
   TBranch *b_global_pot;
   chain->SetBranchAddress("POT_Total",&global_pot,&b_global_pot);
   for(int i =0;i<tot_entries;i++){
     chain->GetEntry(i);
     tot_pot_data+= pot_data;
     global_pot_data+= global_pot;
   
   }
   std::cout<<"Total Data POT: "<<global_pot_data<<std::endl;
   std::cout<<"Used Data POT: "<<tot_pot_data<<std::endl;
   delete chain;
   return tot_pot_data;
   
}

double NukeCCUtilsNSF::setPOTMC(ChainWrapper* _chain){
   TChain *chain = _chain->GetChain();
   double tot_pot_data = 0.0;
   double pot_data = 0.0;
   TBranch *b_pot_data;
   chain->SetBranchAddress("POT_Used",&pot_data,&b_pot_data);
   Int_t tot_entries = chain->GetEntries();
   double global_pot_data = 0.0;
   double global_pot=0.0;
   TBranch *b_global_pot;
   chain->SetBranchAddress("POT_Total",&global_pot,&b_global_pot);
   for(int i =0;i<tot_entries;i++){
     chain->GetEntry(i);
     tot_pot_data+= pot_data;
     global_pot_data+= global_pot;
   
   }
   std::cout<<"Total MC POT: "<<global_pot_data<<std::endl;
   std::cout<<"Used MC POT: "<<tot_pot_data<<std::endl;
   C_global_used_pot_mc = tot_pot_data;
   C_global_tot_pot_mc =  global_pot_data;
   delete chain;
   return tot_pot_data;
   
}

double NukeCCUtilsNSF::GetPOTWeight(){
  ChainWrapper *chain_mc =  GetChainWrapperMCPointer(DefaultCVUniverse::GetPlaylist(),"Meta");
  if(setPOTMC(chain_mc)==0.0)setPOTMC(chain_mc);
  assert(C_global_tot_pot_mc);

  return C_global_used_pot_mc/C_global_tot_pot_mc;
}
/////////////////////////////////////////////Amit Style////////

 NukeCCUtilsNSF::FitBinBoundaries NukeCCUtilsNSF::GetFitBoundaries(int region){
   FitBinBoundaries fit;
   if(region==0){
     fit.region=region;
     fit.y_low = 0.0;
     fit.y_high = 0.2;
     fit.x_low = 1.5;
     fit.x_high = 5.0;
   
   }
   if(region==1){
     fit.region=region;
     fit.y_low = 0.2;
     fit.y_high = 0.4;
     fit.x_low = 1.5;
     fit.x_high = 5.0;    
   
   }
      if(region==2){
     fit.region=region;
     fit.y_low = 0.4;
     fit.y_high = 0.65;
     fit.x_low = 1.5;
     fit.x_high = 5.0;    
   
   }
    if(region==3){
     fit.region=region;
     fit.y_low = 0.65;
     fit.y_high = 0.82;
     fit.x_low = 1.5;
     fit.x_high = 5.0;    
   
   }
   if(region==4){
     fit.region=region;
     fit.y_low = 0.82;
     fit.y_high = 1.0;
     fit.x_low = 1.5;
     fit.x_high = 5.0;    
   
   }
   if(region==5){
     fit.region=region;
     fit.y_low = 1.0;
     fit.y_high = 2.5;
     fit.x_low = 1.5;
     fit.x_high = 5.0;    
   
   }
  if(region==6){
     fit.region=region;
     fit.y_low = 0.0;
     fit.y_high = 0.2;
     fit.x_low = 5.0;
     fit.x_high = 8.0;    
   
   } 
   if(region==7){
     fit.region=region;
     fit.y_low = 0.2;
     fit.y_high = 0.4;
     fit.x_low = 5.0;
     fit.x_high = 8.0;    
   
   }
      if(region==8){
     fit.region=region;
     fit.y_low = 0.4;
     fit.y_high = 0.65;
     fit.x_low = 5.0;
     fit.x_high = 8.0;    
   
   }
    if(region==9){
     fit.region=region;
     fit.y_low = 0.65;
     fit.y_high = 0.82;
     fit.x_low = 5.0;
     fit.x_high = 8.0;    
   
   }
   if(region==10){
     fit.region=region;
     fit.y_low = 0.82;
     fit.y_high = 1.0;
     fit.x_low = 5.0;
     fit.x_high = 8.0;    
   
   }
   if(region==11){
     fit.region=region;
     fit.y_low = 1.0;
     fit.y_high = 2.5;
     fit.x_low = 5.0;
     fit.x_high = 8.0;    
   }
   if(region==12){
     fit.region=region;
     fit.y_low = 0.0;
     fit.y_high = 0.5;
     fit.x_low = 8.0;
     fit.x_high = 15.0;    
   
   }
   
   if(region==13){
     fit.region=region;
     fit.y_low = 0.5;
     fit.y_high = 2.5;
     fit.x_low = 8.0;
     fit.x_high = 15.0;    
   
   } 
   if(region>13 || region<0){
    //std::cout<<"Unphysical fitting region...."<<std::endl;
     fit.region=99;
     fit.y_low = 0.0;
     fit.y_high = 0.0;
     fit.x_low = 0.0;
     fit.x_high = 0.0;   
   
   }     
      
   return fit;
 }
int NukeCCUtilsNSF::GetFitRegion(double pz_val,double pt_val){

 int _region = -1;

 for (int i=0;i<14;i++){
    FitBinBoundaries bound =  GetFitBoundaries(i);
    double pz_low = bound.x_low;
    double pz_high = bound.x_high;
    double pt_low = bound.y_low;
    double pt_high = bound.y_high;
    
    if((pz_val>=pz_low) && (pz_val<pz_high) && (pt_val>=pt_low) && (pt_val<pt_high)){
      _region = i;
      
      break;
    }
 
 }
    return _region;
} 
/////////////////////////////////////////////////////////////// 

void NukeCCUtilsNSF::writePOT( TFile *f ){
  
  double data = 0.0;
  double mc   = 0.0;
  this->getPOT(data,mc);
  TVector2 *pot = new TVector2( data, mc );
   f->WriteTObject( pot, "pot" );
   
 }

//std::map< std::string,std::vector<CVUniverse*>>NukeCCUtilsNSF::GetCVErrorBands(std::map< std::string, std::vector<CVUniverse*> > SystMap){
  // std::map< std::string, std::vector<CVUniverse*> > CVMap;

  // CV
 // CVMap[std::string("CV")] = SystMap[std::string("CV")];
 // return CVMap;
  
//}
  



std::map<std::string,std::vector<CVUniverse*>>NukeCCUtilsNSF::GetErrorBands(ChainWrapper*chain){
  typedef std::map< std::string, std::vector<CVUniverse*> > SystMap;

  SystMap error_bands;

  // CV
  error_bands[std::string("CV")].push_back( new CVUniverse(chain,0) );

  //if RunCodeWithSystematics is not set....dont put these systematics...
  if(RunCodeWithSystematics){
  
  //Flux
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
 
#endif
