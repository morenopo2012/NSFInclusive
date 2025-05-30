//#include "include/CCQENuPlotUtils.h"
//#include "../plot_macros_pub/EventRate/drawUtils.h"
//#include "include/CCQENuBinning.h"
#include "../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
using namespace NUKECC_ANA;
using namespace PlotUtils;

TH2D NormalizeMigrationHistogram( const TH2D* h_migration );

void bookHistos( TFile* file, TH2D** h, string var_name , MnvPlotter *plotter){

  for( unsigned int i = 0; i < nHistos; i++ ){
    h[i] = (TH2D*)file->Get( ( var_name + '_' + names[i] ).c_str() );
    if( h[i] ){ 
      plotter->ApplyAxisStyle( h[i] );
    } else {
      cout << "Couldn't get TH2D from file: " << var_name + '_' + names[i] << endl;
    }
  }

}

int MigrationMatrixPlots( string location)
{
  // read file to get histograms 
  TFile *f = new TFile(Form("%s/Hists_Migration_t3_z26_Nu_minervame1D.root",location.c_str()), "READ" );
  /*if (f->IsZombie() || f->GetListOfKeys()->IsEmpty()){
    Error("MigrationMatrixPlots","Could not get histogram ROOT file or it was empty.");
    return 1;
  }*/
   cout<<"Histname"<<f->GetName()<<endl;
  f->ls();
  MnvPlotter *plotter;
  plotter->SetROOT6Palette(87);
  gStyle->SetNumberContours(200);
   TCanvas* c;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //-------------------------------------------
  // grab event counts before cuts for norm
  //-------------------------------------------
  //TVector2 *evt = (TVector2*)f->Get("n_events");
  //double mc_events = evt->Y();

  //cout<< " Number of MC Events = " <<   mc_events << endl;
  //--------------------------------------------
  // Load histos to plot 
  //--------------------------------------------
  MnvH2D *h_pzmu_ptmu_migration, *h_pzmu_ptmu, *h_ptmu_migration, *h_pzmu_migration,*h_enu_migration[nHistos],*h_q2_migration[nHistos];  
/*  // utils->bookHistos( f, h_pzmu_ptmu_migration, "h_pzmu_ptmu_migration" );
  //utils->bookHistos( f, h_enu_migration, "h_enu_migration" );
  //utils->bookHistos( f, h_q2_migration, "h_q2_migration" );
  bookHistos( f, h_pzmu_migration, "h_pzmu_migration" );
  bookHistos( f, h_ptmu_migration, "h_ptmu_migration" );
*/h_pzmu_ptmu_migration = (MnvH2D*)f->Get("selected_mc_response2d_Emu_Ehad_migration");
  h_pzmu_ptmu = (MnvH2D*)f->Get("h_mc_true_Emu_Ehad");
  h_pzmu_migration = (MnvH2D*)f->Get("selected_Migration_Emu_Ehad");
  h_ptmu_migration = (MnvH2D*)f->Get("selected_Migration_Ehad_Ehad");

  
  //--------------------------------------------
  // Plot Variables
  //--------------------------------------------
  cout << "Plotting Smearing Histograms..." << endl;
  c = new TCanvas("cMig","Smearing");
  
  //------------------------------------------------------
  //Draw entire migration matrix first
  //-------------------------------------------------------
  h_pzmu_ptmu_migration->GetZaxis()->SetTitle("# Normalized Entries");
  h_pzmu_ptmu_migration->GetZaxis()->SetRangeUser(0.,1.); 


  //----------------------------------------------------------
  //Draw back lines to separate PZ bins in each PT space
  //----------------------------------------------------------
  const int nbins_pz = h_pzmu_ptmu->GetNbinsX()+2;
  const int nbins_pt = h_pzmu_ptmu->GetNbinsY()+2;
  //x,y are the migration axis
  const int nbins_x  = h_pzmu_ptmu_migration->GetNbinsX();
  const int nbins_y  = h_pzmu_ptmu_migration->GetNbinsY();
  const double min_x = 0.0; 
  const double min_y = 0.0;
  const double max_x = h_pzmu_ptmu_migration->GetXaxis()->GetBinUpEdge( nbins_x ); 
  const double max_y = h_pzmu_ptmu_migration->GetYaxis()->GetBinUpEdge( nbins_y ); 

  TH2D h_migration_normalized = NormalizeMigrationHistogram( h_pzmu_ptmu_migration );
  h_migration_normalized.GetXaxis()->SetTitle("Reconstructed Emu per Ehad bin");
  h_migration_normalized.GetYaxis()->SetTitle("True Emu  per Ehad  bin");
  h_migration_normalized.GetXaxis()->SetTitleOffset(1.2);
  h_migration_normalized.GetYaxis()->SetTitleOffset(1.2);
  h_migration_normalized.GetZaxis()->SetTitleOffset(1.2);

  h_migration_normalized.Draw("colz");

  cout << nbins_pt << "\t" << nbins_pz << endl;
  
  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(3);
  for(int i=1; i< nbins_pt; ++i ) {
    line.DrawLine( i* nbins_pz, min_y, i*nbins_pz, max_y );
    line.DrawLine( min_x, i* nbins_pz, max_x, i* nbins_pz );
  }
  //  plotter->AddHistoTitle( "Full Migration Matrix", 0.04 );
  c->Print("migration_matrix.png");
  //c->Print("migration_matrix.eps");
  //c->Print("migration_matrix.C");

  //-------------------------------------------------------
  // Draw the PZ migration matrix now 
  //-------------------------------------------------------
  h_pzmu_migration->GetZaxis()->SetTitle("# Normalized Entries");
  h_pzmu_migration->GetZaxis()->SetRangeUser(0.,1.);
  h_pzmu_migration->GetZaxis()->SetTitleOffset(1.2);
  h_pzmu_migration->GetXaxis()->SetTitle("Reconstructed E_{#mu} (GeV)");
  h_pzmu_migration->GetYaxis()->SetTitle("True Muon E_{#mu} (GeV)");

  //-------------------------------------------------------
  // Draw the PZ migration matrix with the corresponding 
  // binwidth representation for each bin 
  //-------------------------------------------------------


  TH2D h_pzmu_migration_normalized = NormalizeMigrationHistogram( h_pzmu_migration );
  h_pzmu_migration_normalized.GetZaxis()->SetRangeUser(0,1.);
  h_pzmu_migration_normalized.Draw("colztext");
  c->Update(); 
  c->Print("Emu_migration_binwidth.png"); 

  //-------------------------------------------------------
  // Draw the PZ migration matrix as a simple matrix 
  // with the smearing number for each bin  
  //-------------------------------------------------------
  cout<<"what is the problem 0"<<endl;
/*  plotter->DrawNormalizedMigrationHistogram( h_pzmu_migration, true, false, false ); 
  c->Update(); 
  c->Print("pzmu_migration_smearing.png"); 
*/
  cout<<"what is the problem 1"<<endl;

  //-------------------------------------------------------
  //Draw the PT migration matrix now 
  //-------------------------------------------------------
  /*h_q2_migration[kCC]->GetZaxis()->SetTitle("# Normalized Entries");
  h_q2_migration[kCC]->GetZaxis()->SetTitleOffset(1.05);
  h_q2_migration[kCC]->GetXaxis()->SetTitle("Reconstructed Q^{2}_{QE} [GeV^{2}]");
  h_q2_migration[kCC]->GetYaxis()->SetTitle("True Q^{2}_{QE} [GeV^{2}]");
  //-------------------------------------------------------
  // Draw the PT migration matrix with the corresponding 
  // binwidth representation for each bin 
  //-------------------------------------------------------
  TH2D h_q2_migration_normalized = NormalizeMigrationHistogram( h_pzmu_migration );
  NukeCC_Binning  *binsDef = new NukeCC_Binning();
  std::vector<double> Q2bins = binsDef->GetDISBins("Emu"); 
  //axis_binning Q2bins = binner->GetQ2BinsGeV();
  Q2bins.bin_edges[0]=0.001;
  TH2D *h_q2_tmp = new TH2D("modified_binning","modified_binning",Q2bins.nbins, &(Q2bins.bin_edges[0]),Q2bins.nbins, &(Q2bins.bin_edges[0]));
  h_q2_tmp->GetZaxis()->SetRangeUser(0,1.0);
  for(int x=0;x<Q2bins.nbins+2;x++){
    for(int y=0;y<Q2bins.nbins+2;y++){
      h_q2_tmp->SetBinContent(x,y,h_q2_migration_normalized.GetBinContent(x,y));
    }
  }
    

  h_q2_tmp->Draw("colz");
  c->SetLogx(1);
  c->SetLogy(1);
  c->Update(); 
  c->Print("q2_migration_binwidth.png"); 
  c->Print("q2_migration_binwidth.eps"); 
  c->Print("q2_migration_binwidth.C"); 
  
  //-------------------------------------------------------
  // Draw the PT migration matrix as a simple matrix  
  // with the smearing number for each bin 
  //-------------------------------------------------------
  plotter->DrawNormalizedMigrationHistogram( h_q2_migration[kCC], true, false, false ); 
  c->Update(); 
  c->Print("q2_migration_smearing.png"); 
  c->Print("q2_migration_smearing.eps"); 
  c->Print("q2_migration_smearing.C"); 
  */
  

  //c->SetLogx(0);
  //c->SetLogy(0);
  //c->Update(); 

  //-------------------------------------------------------
  //Draw the PT migration matrix now 
  //-------------------------------------------------------

  h_ptmu_migration->GetZaxis()->SetTitle("# Normalized Entries");
  h_ptmu_migration->GetZaxis()->SetTitleOffset(1.2);
  h_ptmu_migration->GetXaxis()->SetLabelOffset(0.015);
  h_ptmu_migration->GetYaxis()->SetLabelOffset(0.015);
  h_ptmu_migration->GetXaxis()->SetTitleOffset(1.3);
  h_ptmu_migration->GetYaxis()->SetTitleOffset(1.3);
  h_ptmu_migration->GetXaxis()->SetTitle("Reconstructed E_{recoil} (GeV)");
  h_ptmu_migration->GetYaxis()->SetTitle("True Muon E_{recoil} (GeV)");
  
  //-------------------------------------------------------
  // Draw the PT migration matrix with the corresponding 
  // binwidth representation for each bin 
  //-------------------------------------------------------
  TH2D h_ptmu_migration_normalized = NormalizeMigrationHistogram( h_ptmu_migration );
  h_ptmu_migration_normalized.GetZaxis()->SetRangeUser(0,1.);
  h_ptmu_migration_normalized.Draw("colztext");
  //c->SetLogx();
  //c->SetLogy();
  c->Update(); 
  c->Print("Ehad_migration_binwidth.png"); 


  
  //-------------------------------------------------------
  // Draw the PT migration matrix as a simple matrix  
  // with the smearing number for each bin 
  //-------------------------------------------------------
/*  plotter->DrawNormalizedMigrationHistogram( h_ptmu_migration, true, false, false ); 
  c->Update(); 
  c->Print("ptmu_migration_smearing.png"); 
*/
  //-------------------------------------------------------
  //Draw the ENU migration matrix now 
  //-------------------------------------------------------
  /*  h_enu_migration[kCC]->GetZaxis()->SetTitle("# Normalized Entries");
  h_enu_migration[kCC]->GetZaxis()->SetTitleOffset(1.05);
  h_enu_migration[kCC]->GetXaxis()->SetTitle("Reconstructed E_{#nu,QE} [GeV]");
  h_enu_migration[kCC]->GetYaxis()->SetTitle("True E_{#nu,QE} [GeV]");
  //-------------------------------------------------------
  // Draw the PT migration matrix with the corresponding 
  // binwidth representation for each bin 
  //-------------------------------------------------------
  TH2D h_enu_migration_normalized = NormalizeMigrationHistogram( h_enu_migration[kCC] );
  h_enu_migration_normalized.GetZaxis()->SetRangeUser(0,1.);
  h_enu_migration_normalized.Draw("colz");
  c->Update(); 
  c->Print("enu_migration_binwidth.png"); 
  c->Print("enu_migration_binwidth.eps"); 
  c->Print("enu_migration_binwidth.C"); 
  
  
  //-------------------------------------------------------
  // Draw the PT migration matrix as a simple matrix  
  // with the smearing number for each bin 
  //-------------------------------------------------------
  plotter->DrawNormalizedMigrationHistogram( h_enu_migration[kCC], true, false, false ); 
  c->Update(); 
  c->Print("enu_migration_smearing.png"); 
  c->Print("enu_migration_smearing.eps"); 
  c->Print("enu_migration_smearing.C"); 
*/


  //-------------------------------------------------------------
  //Draw PZ smearings in its corresponding gen. PT_i
  //for the same reco PT_i and the right neighbor (reco PT_(i+1))
  //--------------------------------------------------------------
  gStyle->SetPaintTextFormat("2.1f");
  h_migration_normalized.SetMarkerSize(2.0);
  for(int i=1; i< nbins_pt-1; ++i )
  {
    c->Clear();
    c->Divide(2,1);
    c->cd(1);
    int min = i*nbins_pz+1;
    int max = (i+1)*nbins_pz;
    h_migration_normalized.GetXaxis()->SetRange(min, max);
    h_migration_normalized.GetYaxis()->SetRange(min, max);
    h_migration_normalized.DrawCopy("colz text");
    c->cd(2);
    int min_next = (i+1)*nbins_pz+1;
    int max_next = (i+2)*nbins_pz;
    h_migration_normalized.GetXaxis()->SetRange(min_next, max_next);
    h_migration_normalized.GetYaxis()->SetRange(min, max);
    h_migration_normalized.DrawCopy("colz text");
    c->Print(Form("migration_matrix_pz_genptbin_%i.png",i+1));
    //c->Print(Form("migration_matrix_pz_genptbin_%i.eps",i+1));
    //c->Print(Form("migration_matrix_pz_genptbin_%i.C",i+1));
  }

  f->Close();
  delete f;


  return 0;
};

TH2D NormalizeMigrationHistogram( const TH2D* h_migration )
{

  Int_t nbins = h_migration->GetNbinsX()+2; //under/overflow bins

  //  TMatrixD m_migration(nbins,nbins);
  TH2D tmp(*h_migration);
  tmp.Reset();
  cout << "Normalize with " << nbins << endl;
  //--------------------------------------------------
  // Normalize all x bins for every y bin selected
  //--------------------------------------------------
  for (int y=0; y<nbins+1; ++y) {
    Double_t norm = 0.;
    for (int x=0; x<nbins; ++x) {
      norm += h_migration->GetBinContent(x,y);
    }

    if (fabs(norm)>1E-8) {
      for (int x=0; x<nbins; ++x) {
        double percentage =  1 * h_migration->GetBinContent(x,y) / norm;
	//        m_migration[y][x] = percentage; //yeah that's right  y/x
        tmp.SetBinContent( x, y, percentage);
      }
    }
  }//End of for loop over y bins 
  return tmp;
}

int main( int argc, char *argv[])
{
//  ROOT::Cintex::Cintex::Enable();
  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<< std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t./MigrationMatrixPlots Name_and_path_to_Input_file\n\n"<<
      "\t-Name_and_path_to_Input_file\t =\t Name of and path to the Input file that will be read in" << std::endl; 
    std::cout<<"-----------------------------------------------------------------------------------------------"<< std::endl;
    return 0; 
  }

  // string outdir=argv[1];
  //! Default parameters
  std::vector<std::string> par;
  par.push_back("MigrationMatrixPlots");
  //par.push_back( Form("%s/ana/rootfiles/MigrationHists.root",getenv("MY_NSFNUKECC") ) );
  //par.push_back( Form("%s/ana/rootfiles/MigrationHists.root",getenv("MY_NSFNUKECC") ) );
  par.push_back( Form("../Hists_Migration_t3_z26_Nu_minervame1D.root" ));
  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;
  
  return MigrationMatrixPlots(par[1] );
}
