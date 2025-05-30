#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"
using namespace CCQENU_ANA;


double NodePID(double* nodeEnergy,int n_nodes,CCQENuCuts *cutter, CCQENuEvent* mc, double recoEnergy=0){
  vector<double> means;
  means.push_back(31.302);
  means.push_back(11.418);
  means.push_back(9.769);
  means.push_back(8.675);
  means.push_back(7.949);
  vector<double> sigmas;
  sigmas.push_back(8.997);
  sigmas.push_back(3.075);
  sigmas.push_back(2.554);
  sigmas.push_back(2.484);
  sigmas.push_back(2.232);
  
  double nodeEnergyVal = 0;
  double chi2=0;
  if(n_nodes>5){
    for(int i=0;i<n_nodes;i++){
      if(i==6) break;
      if(i==0) nodeEnergyVal+=nodeEnergy[0];
      else if(i==1) nodeEnergyVal+=nodeEnergy[1];
      else nodeEnergyVal=nodeEnergy[i];
      if(i>=1){
	chi2+=(nodeEnergyVal-means[i-1])*(nodeEnergyVal-means[i-1])/(sigmas[i-1]*sigmas[i-1]);
      }
    }
  }
  
  else{
    bool pass = cutter->passProtonNodeCut(mc);
    if(pass) chi2=0;
    else chi2=75;
  }
  return chi2;
}
int MuonSelectionHists(string filename, string sample, bool makeFluxConstraintHisto, int multiplicity, string playlist, int n_mcfiles = -1, int n_datafiles = -1, int xtype=-1, int ztype=-1)
{
  //---------------------------------------------
  // create file to store histograms
  //---------------------------------------------
  //TFile *f = new TFile( filename.c_str(), "RECREATE" );
  bool typeCheck = xtype==-1 && ztype==-1;//Do normal
  CCQENuUtils  *utils  = new CCQENuUtils( false, makeFluxConstraintHisto );
  utils->setPlaylist(playlist);
  utils->setFluxReweighterPlaylist();
  utils->setAnalysisType(k2DTransverse);
  CCQENuCuts   *cutter = new CCQENuCuts();
  AnaBinning *binner = new AnaBinning();
  CCQENuBinning *minmodbinner = new CCQENuBinning();
  NeutronBlobBinning *neutbinner = new NeutronBlobBinning();
  GeoUtils *geoUtils = new GeoUtils();

  axis_binning thetabins      = minmodbinner->GetMuonThetaBins();
  axis_binning costhetabins   = minmodbinner->GetMuonCosThetaBinsMiniBoone();

  axis_binning muonTbins      = binner->GetMuonEnergyBinsGeV();
  axis_binning muonPbins      = binner->GetMuonEnergyBinsGeV();
  //axis_binning muonPbins      = binner->GetMuonEnergyUniformBinsGeV();
  //axis_binning muonPtbins      = minmodbinner->GetMuonPtUniformBinsGeV();
  axis_binning muonPtbins     = minmodbinner->GetMuonPtBinsGeV();//13
  axis_binning muonPzbins     = minmodbinner->GetMuonPzBinsGeV();//12
  axis_binning Q2bins         = minmodbinner->GetQ2BinsGeV();//19
  axis_binning Enubins        = minmodbinner->GetMuonPzBinsGeV();//12
  axis_binning blobEnergyBinsMeV= minmodbinner->GetBlobEnergyBinsMeV();//1000
  axis_binning dalphatbins    = minmodbinner->GetDalphaT();
  axis_binning dphitbins    = minmodbinner->GetDphiT();
  axis_binning pnbins    = minmodbinner->GetPn();
  axis_binning dptbins    = minmodbinner->GetDeltaPt();
  axis_binning dptxbins    = minmodbinner->GetDeltaPtx();
  axis_binning dptybins    = minmodbinner->GetDeltaPty();
  axis_binning protonThetabins = minmodbinner->GetProtonThetaBins();
  axis_binning protonKineticbins = minmodbinner->GetProtonEnergyBinsGeV();
  axis_binning signbins          = minmodbinner->GetSignedBins();
  axis_binning signeddalphatbins = minmodbinner->GetSignedDeltaAlphaT();
  axis_binning signeddphitbins = minmodbinner->GetSignedDeltaPhiT();


  axis_binning dthetaPerpbins = neutbinner->GetThetaPerpBinsDegree();
  axis_binning dthetaReactbins = neutbinner->GetThetaReactBinsDegree();

  axis_binning dthetaPerpbins_TKI = neutbinner->GetThetaPerpBinsDegree_TKI();
  axis_binning dthetaReactbins_TKI = neutbinner->GetThetaReactBinsDegree_TKI();
  
  axis_binning dalphatbins_TKI    = minmodbinner->GetDalphaT_TKI();
  axis_binning dphitbins_TKI    = minmodbinner->GetDphiT_TKI();
  axis_binning pnbins_TKI    = minmodbinner->GetPn_TKI();
  axis_binning dptbins_TKI    = minmodbinner->GetDeltaPt_TKI();
  axis_binning dptxbins_TKI    = minmodbinner->GetDeltaPtx_TKI();
  axis_binning dptybins_TKI    = minmodbinner->GetDeltaPty_TKI();
  axis_binning protonThetabins_TKI = minmodbinner->GetProtonThetaBins_TKI();
  axis_binning protonKineticbins_TKI = minmodbinner->GetProtonEnergyBinsGeV_TKI();
  axis_binning signeddalphatbins_TKI = minmodbinner->GetSignedDeltaAlphaT_TKI(); 
  axis_binning signeddphitbins_TKI = minmodbinner->GetSignedDeltaPhiT_TKI();




  //---------------------------------------------
  // Book Histograms
  //---------------------------------------------
  //

  //1D

  TH2D *h_dalphatres,*h_dphitres,*h_pnres,*h_dptres,*h_dptxres,*h_dptyres,*h_protonTnres,*h_protonAngleres;
  TH2D *h_dalphatdiff,*h_dphitdiff,*h_pndiff,*h_dptdiff,*h_dptxdiff,*h_dptydiff,*h_protonTndiff,*h_protonAnglediff;
  TH2D *h_dalphatmig,*h_dphitmig,*h_pnmig,*h_dptmig,*h_dptxmig,*h_dptymig,*h_protonTnmig,*h_protonAnglemig;

  h_dalphatres = new TH2D("h_dalphatres","h_dalphatres",360,0,180,200,-2,2);
  h_dphitres = new TH2D("h_dphitres","h_dphitres",360,0,180,200,-2,2);
  h_pnres = new TH2D("h_pnres","h_pnres",400,0,2,200,-2,2);
  h_dptres = new TH2D("h_dptres","h_dptres",400,0,2,200,-2,2);
  h_dptxres = new TH2D("h_dptxres","h_dptxres",800,-2,2,200,-2,2);
  h_dptyres = new TH2D("h_dptyres","h_dptyres",800,-2,2,200,-2,2);
  h_protonTnres = new TH2D("h_protonTnres","h_protonTnres",500,0,1.5,200,-2,2);
  h_protonAngleres = new TH2D("h_protonAngleres","h_protonAngleres",360,0,180,200,-2,2);


  h_dalphatdiff = new TH2D("h_dalphatdiff","h_dalphatdiff",360,0,180,360,-180,180);
  h_dphitdiff = new TH2D("h_dphitdiff","h_dphitdiff",360,0,180,360,-180,180);
  h_pndiff = new TH2D("h_pndiff","h_pndiff",400,0,2,400,-2,2);
  h_dptdiff = new TH2D("h_dptdiff","h_dptdiff",400,0,2,400,-2,2);
  h_dptxdiff = new TH2D("h_dptxdiff","h_dptxdiff",800,-2,2,800,-2,2);
  h_dptydiff = new TH2D("h_dptydiff","h_dptydiff",800,-2,2,800,-2,2);
  h_protonTndiff = new TH2D("h_protonTndiff","h_protonTndiff",500,0,1.5,500,-1.5,1.5);
  h_protonAnglediff = new TH2D("h_protonAnglediff","h_protonAnglediff",360,0,180,360,-180,180);


  h_dalphatmig = new TH2D("h_dalphatmig","h_dalphatmig",360,0,180,360,0,180);
  h_dphitmig = new TH2D("h_dphitmig","h_dphitmig",360,0,180,360,0,180);
  h_pnmig = new TH2D("h_pnmig","h_pnmig",400,0,2,400,0,2);
  h_dptmig = new TH2D("h_dptmig","h_dptmig",400,0,2,400,0,2);
  h_dptxmig = new TH2D("h_dptxmig","h_dptxmig",800,-2,2,800,-2,2);
  h_dptymig = new TH2D("h_dptymig","h_dptymig",800,-2,2,800,-2,2);
  h_protonTnmig = new TH2D("h_protonTnmig","h_protonTnmig",500,0,1.5,500,0,1.5);
  h_protonAnglemig = new TH2D("h_protonAnglemig","h_protonAnglemig",360,0,180,360,0,180);


  
  

  //--------------------------------------------------------------
  // Add Vertical and Lateral Error Bands
  // JO is removing error bands from all 1D because they are all 
  // present in 2D. Plz make projections directly from the 2D. 
  //--------------------------------------------------------------
  TChain *truth_mc = utils->getMCTree("Truth", n_mcfiles );
  
  //---------------------------------------------
  // Get MC Tree
  //---------------------------------------------
  TChain* tree_mc   =  utils->getMCTree("CCQENu", n_mcfiles );
  int entries_mc   = tree_mc->GetEntries();
  utils->setmnvHadronReweightTruthTree(truth_mc);
  utils->setmnvHadronReweightDataTree(tree_mc);
  cout << "MC entries: " << entries_mc << endl;
  CCQENuEvent* mc = new CCQENuEvent( tree_mc );

  //-----------------------------------------------------------
  // Running over event selection criteria for each MC event 
  //-----------------------------------------------------------
  double mc_events = 0., data_events = 0.;
  int mc_nocuts_1track = 0, mc_cut1_1track = 0, mc_cut2_1track =0, mc_cut3_1track = 0, mc_cut4_1track = 0, mc_cut5_1track = 0, mc_cut6_1track = 0, mc_cut6_1track_sideband = 0;
  int mc_nocuts_2track = 0, mc_cut1_2track = 0, mc_cut2_2track =0, mc_cut3_2track = 0, mc_cut4_2track = 0, mc_cut5_2track = 0, mc_cut6_2track = 0, mc_cut7_2track = 0, mc_cut7_2track_sideband = 0;
  int data_nocuts_1track = 0, data_cut1_1track = 0, data_cut2_1track =0, data_cut3_1track = 0, data_cut4_1track = 0, data_cut5_1track = 0, data_cut6_1track = 0, data_cut6_1track_sideband = 0;
  int data_nocuts_2track = 0, data_cut1_2track = 0, data_cut2_2track =0, data_cut3_2track = 0, data_cut4_2track = 0, data_cut5_2track = 0, data_cut6_2track = 0, data_cut7_2track = 0, data_cut7_2track_sideband = 0;
  
  for (int i=0; i<entries_mc ; ++i) {
    if( (i+1)%100000 == 0 ) 
      cout << "Reading MC entry : " << i+1 << " - " << 100*(i+1)/entries_mc << " % " << endl; 
    mc->GetEntry(i);    
    if ( sample == "Signal") {
      if( multiplicity==0) { //1+2tracks
	if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passCCQESelection(mc) )             continue;
      }
      else if( multiplicity==1) { //1 track only 
	if( !cutter->passInteractionVertex( mc ) )       continue;	mc_nocuts_1track++; 
        if( !cutter->passDeadTimeCuts( mc ) )            continue;	mc_cut1_1track++; 
        if( !cutter->passNuHelicityCut( mc ) )   	 continue;  
        if( !(mc->multiplicity==1) )                     continue;	mc_cut3_1track++; 
        if( !cutter->passImprovedMichelCut( mc ) )       continue;	mc_cut4_1track++; 
	if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( !cutter->passNBlobs( mc ) )                  continue;	mc_cut5_1track++; 
	if ( !cutter->passTrackAngleCut( mc ) )          continue;	mc_cut6_1track++; 
      }
      else if( multiplicity==2) { //2 tracks only 
	if( !cutter->passInteractionVertex( mc ) )       continue; 	mc_nocuts_2track++;  
	if( !cutter->passDeadTimeCuts( mc ) )            continue;	mc_cut1_2track++; 
        if( !cutter->passNuHelicityCut( mc ) )           continue;	mc_cut2_2track++; 
        if( !(mc->multiplicity==2) )                     continue;	mc_cut3_2track++; 
        if( !cutter->passImprovedMichelCut( mc ) )       continue;	mc_cut5_2track++; 
	if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( !cutter->passNBlobs( mc ) )                  continue;	mc_cut6_2track++; 
	if( !cutter->passTrackAngleCut( mc ) )           continue;	mc_cut7_2track++; 
	if( !cutter->passProtonContainment(mc,1100,false)) continue;
      }
    }
    else if( sample == "MichelSideBand") {
      if( multiplicity==0) {//1+2tracks
	if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !cutter->passMultiplicity( mc ) )            continue;
	if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( !cutter->passNBlobs( mc ) )                  continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
	if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
      else if( multiplicity==1) { //1 track only 
	if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==1) )                     continue;
	if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( !cutter->passNBlobs( mc ) )                  continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
	if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
      else if( multiplicity==2) { //2 tracks only 
	if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==2) )                     continue;
	if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( !cutter->passNBlobs( mc ) )                  continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
	if( !cutter->passTrackAngleCut( mc ) )           continue;
	if( !cutter->passProtonContainment(mc,1100,false)) continue;
      }
    }
    else if( sample == "BlobSideBand") {
      if( multiplicity==0) {//1+2tracks
	if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !cutter->passMultiplicity( mc ) )            continue;
	if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( !cutter->passImprovedMichelCut( mc ) )       continue;
	if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
      else if( multiplicity==1) { //1 track only 
	if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==1) )                     continue;
	if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( !cutter->passImprovedMichelCut( mc ) )       continue;
	if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
      else if( multiplicity==2) { //2 tracks only 
	if( !cutter->passInteractionVertex( mc ) )       continue;	
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==2) )                     continue;
	if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( !cutter->passImprovedMichelCut( mc ) )       continue;
	if( !cutter->passTrackAngleCut( mc ) )           continue;
	if( !cutter->passProtonContainment(mc,1100,false)) continue;
      }
    }
    else if( sample == "MicBlobSideBand") {
      if( multiplicity==0) {//1+2tracks
        if( !cutter->passInteractionVertex( mc ) )       continue;
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !cutter->passMultiplicity( mc ) )            continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
        if( !cutter->passTrackAngleCut( mc ) )           continue;
      }
      else if( multiplicity==1) { //1 track only 
        if( !cutter->passInteractionVertex( mc ) )       continue;
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==1) )                     continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
        if( !cutter->passTrackAngleCut( mc ) )          continue;
      }
      else if( multiplicity==2) { //2 tracks only      
        if( !cutter->passInteractionVertex( mc ) )       continue;
        if( !cutter->passDeadTimeCuts( mc ) )            continue;
        if( !cutter->passNuHelicityCut( mc ) )           continue;
        if( !(mc->multiplicity==2) )                     continue;
        if( !cutter->passExtraTracksProtonsCut ( mc ) )  continue;
        if( cutter->passNBlobs( mc ) )                   continue;
        if( cutter->passImprovedMichelCut( mc ) )        continue;
        if( !cutter->passTrackAngleCut( mc ) )           continue;
	if( !cutter->passProtonContainment(mc,1100,false)) continue;
      }
    } else {
      cout<<"Wrong Sample Selection."<<endl;
      cout<<"Select: Signal, MichelSideBand, BlobSideBand, MicBlobSideBand"<<endl;
      continue;
    }
    if(!cutter->passTrueCCQELike(mc)) continue;
    double protonAngle = mc->CCQENu_proton_theta*180./3.14159;
    if(protonAngle>70) continue;
    double trackchi2 = NodePID(mc->CCQENu_proton_nodes_nodesNormE,mc->CCQENu_proton_nodes_nodesNormE_sz,cutter,mc,0);
    if(trackchi2>10) continue;
      
    //Get CVweight
    double wgt = utils->GetCVWeight( mc , sample);    
    mc_events+=wgt;

    ///////////////////////////////////
    // RECO SECTION//
    ///////////////////////////////////


    //-----------
    //Muon Variables
    double muon_theta     = mc->muon_theta_allNodes[20]; // get theta from 20th node, note NX kludge fixes this...
    double cos_muon_theta =  cos(muon_theta);
    double muon_T         =  mc->CCQENu_muon_T / pow(10,3); //GeV
    double reco_muon_px   =  mc->CCQENu_leptonE[0];
    double reco_muon_py   =  mc->CCQENu_leptonE[1];
    double reco_muon_pz   =  mc->CCQENu_leptonE[2];
    double muon_p         = utils->GetTotalMomentum( reco_muon_px, reco_muon_py, reco_muon_pz ) / pow(10,3); //GeV
    //New method (maybe temporary)
    double muon_pt_beam = sin(muon_theta)*muon_p;
    double muon_pz_beam = cos(muon_theta)*muon_p;

    //Proton Variables
    TVector3 protonVect(mc->CCQENu_proton_Px_fromdEdx,mc->CCQENu_proton_Py_fromdEdx,mc->CCQENu_proton_Pz_fromdEdx);//MeV
    double protonMom = protonVect.Mag();//MeV
    double protonTn = TMath::Sqrt(protonMom*protonMom+utils->M_p*utils->M_p)-utils->M_p;//Done in MeV

    
    //Convert to degrees, GeV
    muon_theta*= 180. / 3.14159;
    reco_muon_px /= pow(10,3);
    reco_muon_py /= pow(10,3);
    reco_muon_pz /= pow(10,3);
    protonVect *= pow(10,-3);
    protonTn /= pow(10,3);


    double reco_dalphat = utils->GetDeltaAlphaT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_dphit = utils->GetDeltaPhiT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_pn = utils->GetPn(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_dpt = utils->GetDeltaPt(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_dptx = utils->GetDeltaPtx(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_dpty = utils->GetDeltaPty(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_sign = utils->GetMuonProtonSign(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_signed_dalphat = utils->GetSignedDeltaAlphaT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
    double reco_signed_dphit = utils->GetSignedDeltaPhiT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);

    const double bindingE = 34.0;
    /////////////////////////////
    //True Section
    ////////////////////////////

    //Reco object
    TVector3 true_proton_p_comp_reco(mc->proton_prong_4p[1],mc->proton_prong_4p[2],mc->proton_prong_4p[3]);
    



    TVector3 true_proton_p_comp = utils->GetHighestTrueProtonMomentumComp(mc);
    double true_proton_px = true_proton_p_comp.X();
    double true_proton_py = true_proton_p_comp.Y();
    double true_proton_pz = true_proton_p_comp.Z();




    double true_proton_P = TMath::Sqrt(true_proton_px*true_proton_px+
				       true_proton_py*true_proton_py+
				       true_proton_pz*true_proton_pz);

    double true_protonTn = TMath::Sqrt(true_proton_P*true_proton_P+(utils->M_p)*(utils->M_p))-utils->M_p;
    double true_protonAngle = utils->GetTheta(true_proton_px,true_proton_py,true_proton_pz)*180.0/TMath::Pi();


    if(true_protonTn>600) continue; // Cut for the inelastic ones which just die off at a super high rate and get reco'd at halfish the energy.

    double true_muon_px   = mc->mc_primFSLepton[0];
    double true_muon_py   = mc->mc_primFSLepton[1];
    double true_muon_pz   = mc->mc_primFSLepton[2];
    double true_muon_E    = mc->mc_primFSLepton[3];
    double true_muon_p    = utils->GetTotalMomentum( true_muon_px, true_muon_py, true_muon_pz ); //GeV
    double true_theta     = utils->GetTheta( true_muon_px, true_muon_py, true_muon_pz )*180.0/TMath::Pi();
    double true_muon_pt_beam   = utils->GetTransverseMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV 
    double true_muon_pz_beam   = utils->GetLongitudinalMomentumWRTBeam( true_muon_px, true_muon_py, true_muon_pz ); //GeV
  

    //Scales
    true_proton_px /= pow(10,3);
    true_proton_py /= pow(10,3);
    true_proton_pz /= pow(10,3);
    true_proton_P  /= pow(10,3);
    true_protonTn        /= pow(10,3);
    true_muon_px   /= pow(10,3);
    true_muon_py   /= pow(10,3);
    true_muon_pz   /= pow(10,3);
    true_muon_E    /= pow(10,3);
    true_muon_p    /= pow(10,3);
    true_muon_pt_beam /= pow(10,3);
    true_muon_pz_beam /= pow(10,3);
    true_proton_p_comp_reco*= pow(10,-3);




    double true_dalphat   = utils->GetDeltaAlphaT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dphit     = utils->GetDeltaPhiT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_pn        = utils->GetPn(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dpt       = utils->GetDeltaPt(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dptx      = utils->GetDeltaPtx(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_dpty      = utils->GetDeltaPty(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_sign      = utils->GetMuonProtonSign(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_signed_dalphat = utils->GetSignedDeltaAlphaT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
    double true_signed_dphit   = utils->GetSignedDeltaPhiT(true_proton_px,true_proton_py,true_proton_pz,true_muon_px,true_muon_py,true_muon_pz);
 

    if(true_proton_P>1.2 || true_proton_P<0.45 || true_protonAngle>70) continue;//HE proton determined 


    double cheat_proton_dalphat = utils->GetDeltaAlphaT(true_proton_px,true_proton_py,true_proton_pz,reco_muon_px,reco_muon_py,reco_muon_pz);
    double cheat_muon_dalphat   = utils->GetDeltaAlphaT(protonVect.X(),protonVect.Y(),protonVect.Z(),true_muon_px,true_muon_py,true_muon_pz);

    if(true_dalphat<100 && reco_dalphat>150){
      cout << true_proton_p_comp_reco.X() << "\t" << true_proton_p_comp_reco.Y() << "\t" << true_proton_p_comp_reco.Z() << "\t0\t0\t" 
	   << true_proton_px << "\t" <<true_proton_py << "\t" << true_proton_pz << "\t0\t0\t" 
	   <<  protonVect.X() << "\t" << protonVect.Y() << "\t" << protonVect.Z() << "\t0\t0\t" 
	   << true_muon_px << "\t" << true_muon_py << "\t" << true_muon_pz << "\t0\t0\t"  
	   << reco_muon_px << "\t" << reco_muon_py << "\t" << reco_muon_pz << "\t0\t0\t"  
	   << true_dalphat << "\t" << reco_dalphat << "\t" << cheat_proton_dalphat << "\t" << cheat_muon_dalphat << endl;    
    }
    

    //-------------------------------------------------------------------
    // Check that the CV of events pass the proton and recoil cuts 
    // This is because we have systematics due to these cuts 
    // and a change in their values can affect the composition 
    // of the selected sample. 
    //-------------------------------------------------------------------
    bool doesCVpass_SingleProtonCut = (mc->multiplicity==1)? true : cutter->passSingleProtonCut( mc, 0, 0 );
    bool doesCVpass_ExtraProtonsCut = (mc->multiplicity==1)? true : cutter->passExtraProtonsCut( mc, NULL, 0 );
    bool doesCVpass_RecoilCut = cutter->passSignalFunction( mc, 0, 0 ); 

    if( doesCVpass_SingleProtonCut && doesCVpass_ExtraProtonsCut && doesCVpass_RecoilCut ) {     


      h_dalphatres->Fill(true_dalphat,1-reco_dalphat/true_dalphat,wgt);
      h_dphitres->Fill(true_dphit,1-reco_dphit/true_dphit,wgt);
      h_pnres->Fill(true_pn,1-reco_pn/true_pn,wgt);
      h_dptres->Fill(true_dpt,1-reco_dpt/true_dpt,wgt);
      h_dptxres->Fill(true_dptx,1-reco_dptx/true_dptx,wgt);
      h_dptyres->Fill(true_dpty,1-reco_dpty/true_dpty,wgt);
      h_protonTnres->Fill(true_protonTn,1-protonTn/true_protonTn,wgt);
      h_protonAngleres->Fill(true_protonAngle,1-protonAngle/true_protonAngle,wgt);

      h_dalphatdiff->Fill(true_dalphat,reco_dalphat-true_dalphat,wgt);
      h_dphitdiff->Fill(true_dphit,reco_dphit-true_dphit,wgt);
      h_pndiff->Fill(true_pn,reco_pn-true_pn,wgt);
      h_dptdiff->Fill(true_dpt,reco_dpt-true_dpt,wgt);
      h_dptxdiff->Fill(true_dptx,reco_dptx-true_dptx,wgt);
      h_dptydiff->Fill(true_dpty,reco_dpty-true_dpty,wgt);
      h_protonTndiff->Fill(true_protonTn,protonTn-true_protonTn,wgt);
      h_protonAnglediff->Fill(true_protonAngle,protonAngle-true_protonAngle,wgt);

      h_dalphatmig->Fill(true_dalphat,reco_dalphat,wgt);
      h_dphitmig->Fill(true_dphit,reco_dphit,wgt);
      h_pnmig->Fill(true_pn,reco_pn,wgt);
      h_dptmig->Fill(true_dpt,reco_dpt,wgt);
      h_dptxmig->Fill(true_dptx,reco_dptx,wgt);
      h_dptymig->Fill(true_dpty,reco_dpty,wgt);
      h_protonTnmig->Fill(true_protonTn,protonTn,wgt);
      h_protonAnglemig->Fill(true_protonAngle,protonAngle,wgt);


    }
  } //End of for loop over MC entries 
  
  cout << "Finished looping over all MC entries. Number of weighted MC events satisfying all selection criteria = " << mc_events << endl; 
  delete mc; 
  delete tree_mc; 
  
  //---------------------------------------------
  // Get Data Tree
  //---------------------------------------------
  if (n_datafiles==0) {
    cout<<"Data files won't be included"<<endl;
  } else {
    TChain* tree_data  = utils->getDataTree("CCQENu", n_datafiles );
    int entries_data   = tree_data->GetEntries();

    cout << "Data entries: " << entries_data << endl;
    bool isData = true;
    CCQENuEvent* data = new CCQENuEvent( tree_data, isData );

  //---------------------------------------------
  // Fill DATA Histograms
  //---------------------------------------------
    for (int i=0; i<entries_data ; ++i) {
      if( (i+1)%10000 == 0 ) 
	cout << "Reading DATA entry : " << i+1 << " - " << 100*(i+1)/entries_data << " % " << endl; 
      break;
      data->GetEntry(i);

      if( sample == "Signal") {
	if( multiplicity==0) { //1+2tracks
	  if( !cutter->passInteractionVertex( data ) )     continue;	
	  if( !cutter->passCCQESelection( data ) )         continue;
	}
	else if( multiplicity==1) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	  data_nocuts_1track++; 
	  if( !cutter->passDeadTimeCuts( data ) )           continue;	  data_cut1_1track++; 
	  if( !cutter->passNuHelicityCut( data ) ) 	    continue;     data_cut2_1track++; 
	  if( !(data->multiplicity==1) )                    continue;	  data_cut3_1track++; 
	  if( !cutter->passImprovedMichelCut( data ) )      continue; 	  data_cut4_1track++; 
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( !cutter->passNBlobs( data ) )                 continue;	  data_cut5_1track++; 
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;	  data_cut6_1track++; 
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==2) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	  data_nocuts_2track++; 
	  if( !cutter->passDeadTimeCuts( data ) )           continue;	  data_cut1_2track++; 
	  if( !cutter->passNuHelicityCut( data ) )          continue;	  data_cut2_2track++; 
	  if( !(data->multiplicity==2) )                    continue;	  data_cut3_2track++; 
	  if( !cutter->passSingleProtonCut( data ) )        continue;
	  if( !cutter->passExtraProtonsCut( data ) )        continue;	  data_cut4_2track++; 
	  if( !cutter->passImprovedMichelCut( data ) )      continue;	  data_cut5_2track++; 
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( !cutter->passNBlobs( data ) )                 continue;	  data_cut6_2track++; 
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;	  data_cut7_2track++; 
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	  if( !cutter->passProtonContainment(data,1100,false)) continue;
	  if( !cutter->passProtonNodeCut(data) )             continue;

	}
      }
      else if( sample == "MichelSideBand") {
	if( multiplicity==0) {//1+2tracks
	  bool passProtonCut = (data->multiplicity==1)? true : cutter->passSingleProtonCut( data );
	  bool pass_ExtraProtonsCut = (data->multiplicity==1)? true : cutter->passExtraProtonsCut( data );
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !cutter->passMultiplicity( data ) )           continue;
	  if( !passProtonCut )                              continue;
	  if( !pass_ExtraProtonsCut )                       continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( !cutter->passNBlobs( data ) )                 continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==1) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==1) )                    continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( !cutter->passNBlobs( data ) )                 continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==2) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==2) )                    continue;
	  if( !cutter->passSingleProtonCut( data ) )        continue;
	  if( !cutter->passExtraProtonsCut( data ) )        continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( !cutter->passNBlobs( data ) )                 continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if ( !cutter->passTrackAngleCut( data ) )         continue;
	  if( !cutter->passProtonContainment(data,1100,false)) continue;
	  if( !cutter->passProtonNodeCut(data) )             continue;
	}
      }

      else if( sample == "BlobSideBand") {
	if( multiplicity==0) {//1+2tracks
	  bool passProtonCut = (data->multiplicity==1)? true : cutter->passSingleProtonCut( data );
	  bool pass_ExtraProtonsCut = (data->multiplicity==1)? true : cutter->passExtraProtonsCut( data );
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !cutter->passMultiplicity( data ) )           continue;
	  if( !passProtonCut )                              continue;
	  if( !pass_ExtraProtonsCut )                       continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( !cutter->passImprovedMichelCut( data ) )      continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==1) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==1) )                    continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( !cutter->passImprovedMichelCut( data ) )      continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if ( !cutter->passTrackAngleCut( data ) )         continue;
	}
	else if( multiplicity==2) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==2) )                    continue;
	  if( !cutter->passSingleProtonCut( data ) )        continue;
	  if( !cutter->passExtraProtonsCut( data ) )        continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( !cutter->passImprovedMichelCut( data ) )      continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	  if( !cutter->passProtonContainment(data,1100,false)) continue;
	  if( !cutter->passProtonNodeCut(data) )             continue;
	}
      }
      else if( sample == "MicBlobSideBand") {
	if( multiplicity==0) {//1+2tracks
	  bool passProtonCut = (data->multiplicity==1)? true : cutter->passSingleProtonCut( data );
	  bool pass_ExtraProtonsCut = (data->multiplicity==1)? true : cutter->passExtraProtonsCut( data );
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !cutter->passMultiplicity( data ) )           continue;
	  if( !passProtonCut )                              continue;
	  if( !pass_ExtraProtonsCut )                       continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==1) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==1) )                    continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	}
	else if( multiplicity==2) {
	  if( !cutter->passInteractionVertex( data ) )      continue;	
	  if( !cutter->passDeadTimeCuts( data ) )           continue;
	  if( !cutter->passNuHelicityCut( data ) )          continue;
	  if( !(data->multiplicity==2) )                    continue;
	  if( !cutter->passSingleProtonCut( data ) )        continue;
	  if( !cutter->passExtraProtonsCut( data ) )        continue;
	  if( !cutter->passExtraTracksProtonsCut ( data ) ) continue;
	  if( cutter->passNBlobs( data ) )                  continue;
	  if( cutter->passImprovedMichelCut( data ) )       continue;
	  if( !cutter->passSignalFunction( data, 0, 0 ) )   continue;
	  if( !cutter->passTrackAngleCut( data ) )          continue;
	  if( !cutter->passProtonContainment(data,1100,false)) continue;
	  if( !cutter->passProtonNodeCut(data) )             continue;
	}
      }
      else{
	cout << "You gave an invalid sample" << endl;
	return 1;
      }

      data_events++; 
    //-----------
    // This is the old method. Because of angular bias introduced from hadronic energy we need to pick out the right angle from a downstream node
    // New method until general reco is fixed is to get the reco P, a downstream theta (already corrected to the beam direction), and calculate
    // the pt and pz from these two quantities
    //-----------
      
      //Muon Variables
      double muon_theta = data->muon_theta_allNodes[20];// get theta from 20th node
      double cos_muon_theta =  cos(muon_theta);
      double muon_T = data->CCQENu_muon_T / pow(10,3); //GeV
      double reco_muon_px   =  data->CCQENu_leptonE[0];
      double reco_muon_py   =  data->CCQENu_leptonE[1];
      double reco_muon_pz   =  data->CCQENu_leptonE[2];
      double muon_p    = utils->GetTotalMomentum( reco_muon_px, reco_muon_py, reco_muon_pz) / pow(10,3); //GeV
      double muon_pt_beam = sin(muon_theta)*muon_p;
      double muon_pz_beam = cos(muon_theta)*muon_p;

      //Proton Variables
      TVector3 protonVect(data->CCQENu_proton_Px_fromdEdx,data->CCQENu_proton_Py_fromdEdx,data->CCQENu_proton_Pz_fromdEdx);//MeV
      double protonAngle = data->CCQENu_proton_theta*180./3.14159;
      double protonMom = protonVect.Mag();//MeV
      double protonTn = TMath::Sqrt(protonMom*protonMom+utils->M_p*utils->M_p)-utils->M_p;//MeV


      if(protonAngle>70) continue;
      double trackchi2 = NodePID(data->CCQENu_proton_nodes_nodesNormE,data->CCQENu_proton_nodes_nodesNormE_sz,cutter,data,0);
      if(trackchi2>10) continue;


      //Convert to degrees, GeV
      muon_theta*= 180. / 3.14159;
      reco_muon_px /= pow(10,3);
      reco_muon_py /= pow(10,3);
      reco_muon_pz /= pow(10,3);
      protonVect *= pow(10,-3);
      protonTn /= pow(10,3);

      //Transverse variables all inputs are GeV
      double reco_dalphat = utils->GetDeltaAlphaT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_dphit = utils->GetDeltaPhiT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_pn = utils->GetPn(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_dpt = utils->GetDeltaPt(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_dptx = utils->GetDeltaPtx(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_dpty = utils->GetDeltaPty(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_sign = utils->GetMuonProtonSign(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_signed_dalphat = utils->GetSignedDeltaAlphaT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);
      double reco_signed_dphit = utils->GetSignedDeltaPhiT(protonVect.X(),protonVect.Y(),protonVect.Z(),reco_muon_px,reco_muon_py,reco_muon_pz);

      //Directional "transverse" variables, Tejin
      // beam_bias default to 0, using this method allows calculation of angle uncertainties. 
      // 1. Define the beam direction, and muon 4P
      double beam_bias = 0;
      XYZVector beam = geoUtils->BeamAxis( beam_bias ); 
      XYZTVector reco_muon_4P( reco_muon_px, reco_muon_py, reco_muon_pz ,data->CCQENu_muon_E/1e3);
      
      // 2. Compute expected stationary nucleon scattering momentum
      // XYZTVector ComputeExpectedNucleon ( XYZVector &nu, XYZTVector &muon, double ISMass = 0.93827231,double FSMass =   0.93956563, double BindingE = 0.00 ); // vbar p --> mu n, Initial Mass, Final state mass
      double binding_e = 0; // theoretically we should include binding E for carbon.. but I'm comparing to hydrogen in anti-nu, so probably set it to 0.
      double is_neutron_mass = 0.9395654133;
      double fs_proton_mass = 0.93827231;

      XYZTVector expectedProton4P = geoUtils->ComputeExpectedNucleon( beam, reco_muon_4P, is_neutron_mass, fs_proton_mass, binding_e);
      XYZVector expectedProton3P = (XYZVector) expectedProton4P.Vect();

      // 3. 
      // ComputeNeutronAngularVars( XYZVector &nu, XYZVector &expVec, XYZVector &targetVec ), variables are : dThetaP, dThetaR, dThetaTotal,  dPp, dPr, inferred dPp and dPr
      XYZVector protonVect_XYZ( protonVect.X(), protonVect.Y(), protonVect.Z() );
      vector<double> reco_vars = geoUtils->ComputeNeutronAngularVars( beam, expectedProton3P, protonVect_XYZ );
      double reco_dThetaP = reco_vars[0]*180/TMath::Pi();
      double reco_dThetaR = reco_vars[1]*180/TMath::Pi();
      double reco_dTheta = reco_vars[2];
      double reco_dPp = reco_vars[3];
      double reco_dPr = reco_vars[4];
      double reco_dPp_infer = reco_vars[5];
      double reco_dPr_infer = reco_vars[6];

      //I feel there is also a need  to include the left-right uncertainty of the beam, since LR-direction is pretty important....




      map<int,double> varmap;
      varmap[0] = muon_pz_beam;//muonPzbins;
      varmap[1] = reco_dalphat;//dalphatbins;
      varmap[2] = reco_dphit;//dphitbins;
      varmap[3] = reco_pn;//pnbins;
      varmap[4] = reco_dpt;//dptbins;
      varmap[5] = reco_dptx;//dptxbins;
      varmap[6] = reco_dpty;//dptybins;
      varmap[7] = protonTn;//protonKineticbins;
      varmap[8] = protonAngle;//protonThetabins;
      varmap[9] = reco_sign;//signbins;
      varmap[10] = reco_signed_dalphat;//signeddalphatbins;
      varmap[11] = reco_signed_dphit;//signeddphitbins;
      varmap[12] = reco_dThetaP;  
      varmap[13] = reco_dThetaR;  
      varmap[14] = reco_dTheta;
      varmap[15] = reco_dPp;      
      varmap[16] = reco_dPr;     
      varmap[17] = reco_dPp_infer;
      varmap[18] = reco_dPr_infer;


      
    } //End of for loop over DATA entries 
    
    cout << "Finished looping over all DATA entries. Number of DATA events satisfying all selection criteria = " << data_events << endl; 
    delete data; 
    delete tree_data; 
  }

  //Print out some information
  cout<<" --------------------------------------------- " << endl;
  cout<<" MC " << endl;
  cout<<" --------------------------------------------- " << endl;
  cout<<" ****************** " << endl;
  cout<<" MC 1 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut6_sideband = %i", mc_nocuts_1track, mc_cut1_1track, mc_cut2_1track, mc_cut3_1track, mc_cut4_1track, mc_cut5_1track, mc_cut6_1track, mc_cut6_1track_sideband) << endl;
  cout<<" ****************** " << endl;
  cout<<" MC 2 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut7 = %i, cut7_sideband = %i", mc_nocuts_2track, mc_cut1_2track, mc_cut2_2track, mc_cut3_2track, mc_cut4_2track, mc_cut5_2track, mc_cut6_2track, mc_cut7_2track, mc_cut7_2track_sideband) << endl;
  cout<<" --------------------------------------------- " << endl;
  cout<<" DATA " << endl;
  cout<<" --------------------------------------------- " << endl;
  cout<<" ****************** " << endl;
  cout<<" DATA 1 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut6_sideband = %i", data_nocuts_1track, data_cut1_1track, data_cut2_1track, data_cut3_1track, data_cut4_1track, data_cut5_1track, data_cut6_1track, data_cut6_1track_sideband) << endl;
  cout<<" ****************** " << endl;
  cout<<" DATA 2 Track Events " << endl;
  cout<<" ****************** " << endl;
  cout<< Form("nocuts = %i; cut1 = %i, cut2 = %i, cut3 = %i, cut4 = %i, cut5 = %i, cut6 = %i, cut7 = %i, cut7_sideband = %i", data_nocuts_2track, data_cut1_2track, data_cut2_2track, data_cut3_2track, data_cut4_2track, data_cut5_2track, data_cut6_2track, data_cut7_2track, data_cut7_2track_sideband) << endl;
  
  //==================================================================
  // Create ROOT file to store histograms
  // Write to file all the created histograms 
  //==================================================================
  //  if(!typeCheck) filename = filename+"_"+axis_name[xtype]+"_"+axis_name[ztype]+".root";
  TFile *f = new TFile( filename.c_str(), "RECREATE" );
   
  //Write MC and DATA POT to file
  utils->writePOT( f );
  
  //Write multiplicity and sample definitions
  TMap *sampleMap = new TMap(2,0);
  string multiplicity_label = (multiplicity==0)? "alltracks": Form("%i_track",multiplicity);
  sampleMap->Add(new TObjString("multiplicity"), new TObjString(multiplicity_label.c_str()) );
  sampleMap->Add(new TObjString("sample"), new TObjString(sample.c_str()) );
  f->WriteTObject( sampleMap, "sample");
  
  //Write DATA and MC number of events
  TVector2 *evt = new TVector2( data_events, mc_events );
  f->WriteTObject( evt, "n_events" );

  //Write out whether this ROOT file contains histograms that will be used for the flux constraint  
  TParameter<bool> *fluxConstrain = new TParameter<bool>( "fluxConstraintHisto", makeFluxConstraintHisto ); 
  f->WriteTObject( fluxConstrain ); 
  
  f->cd();
  

  h_dalphatres->Write();
  h_dphitres->Write();
  h_pnres->Write();
  h_dptres->Write();
  h_dptxres->Write();
  h_dptyres->Write();
  h_protonTnres->Write();
  h_protonAngleres->Write();

  h_dalphatdiff->Write();
  h_dphitdiff->Write();
  h_pndiff->Write();
  h_dptdiff->Write();
  h_dptxdiff->Write();
  h_dptydiff->Write();
  h_protonTndiff->Write();
  h_protonAnglediff->Write();


  h_dalphatmig->Write();
  h_dphitmig->Write();
  h_pnmig->Write();
  h_dptmig->Write();
  h_dptxmig->Write();
  h_dptymig->Write();
  h_protonTnmig->Write();
  h_protonAnglemig->Write();
  
  f->Close();
  delete f;
  delete utils;
  delete cutter; 
  delete binner;
  delete minmodbinner; 

  return 0;
}

int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./MuonSelectionHists  Name_and_path_to_Output_file  Selection_Signal  Make_Flux_Constraint_Histo  Multiplicity  Number_MC_files  Number_DATA_files \n\n"<<
      "\t-Name_and_path_to_Output_file\t =\t Name of and path to the Output ROOT file that will be created \n"<<
      "\t-Playlist\t =\t Name of the playlist you want to run over e.g. minerva1 \n"<<
      "\t-Selection_Signal\t =\t Can be: Signal, MichelSideBand, SideBand (this is:non-vtx vs Q2 sideband) \n"<<
      "\t-Make_Flux_Constraint_Histo\t =\t If TRUE Enter 1; If FALSE Enter 0 \n"<<
      "\t-Multiplicity\t =\t Enter 0, 1 or 2. Here 0: 1 or 2 tracks; 1: 1-track events only; 2: 2-track events only \n"<<
      "\t-Number_MC_files\t =\t Number of MonteCarlo files. To use all files, set this to: -1 \n"<<
      "\t-Number_DATA_files\t =\t Number of Data files. To use all files, set this to: -1\n" <<
      "\t-Variable combo\t = \tChoose TKI vs TKI combinations, default is off (-1)" << std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }
  
  //! Default parameters
  std::vector<std::string> par;
  par.push_back("MuonSelectionHists");
  par.push_back( Form("%s/ana/rootfiles/MuonSelectionHists.root",getenv("CCQENUROOT") ) );
  par.push_back("minerva1");
  par.push_back("Signal");
  par.push_back("0");
  par.push_back("0");
  par.push_back("-1");
  par.push_back("-1");
  par.push_back("-1");//no special variable combinations
  par.push_back("-1");

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  bool fluxConstraintHistoNeeded = ( par.at(4) == "1" ) ? true : false; 

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

  return MuonSelectionHists(par[1], par[3], fluxConstraintHistoNeeded, atoi(par[5].c_str()), par[2],atoi(par[6].c_str()), atoi(par[7].c_str()), atoi(par[8].c_str()), atoi(par[9].c_str()) );
}
