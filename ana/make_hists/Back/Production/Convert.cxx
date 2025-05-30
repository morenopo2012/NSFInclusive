#include "TFile.h"
#include "TMath.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "../../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
//#include "syst_common.h"
#include "TArrayD.h"
#include <iostream>
                              
#include "Math/Factory.h"     
#include "Math/Functor.h"     
#include "Math/Minimizer.h"   

#include "RooUnfold/RooUnfoldBayes.h"
#include "RooUnfold/RooUnfoldSvd.h"
#include "RooUnfold/RooUnfoldTUnfold.h"
#include "RooUnfold/RooUnfoldInvert.h"
#include "RooUnfold/RooUnfoldBinByBin.h"
#include "RooUnfold/RooUnfoldResponse.h"
#include "RooUnfold/RooUnfold.h"
#include "MinervaUnfold/MnvResponse.h"
#include "MinervaUnfold/MnvUnfold.h"

#include "PlotUtils/MacroUtil.h" 
//#include "include/CVUniverse_faiza.h"
#include "../../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "../../include/Variable.h"
//#include "../../include/Variable_faiza.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "../../../NUKECCSRC/ana_common/include/LateralSystematics.h"
#include <iostream>
#include <stdlib.h>
//#include "Cintex/Cintex.h"
#include "TParameter.h"
#include "../../include/plotting_functions.h"
#ifndef __CINT__
#endif

using namespace NUKECC_ANA;
using namespace PlotUtils;
using namespace std;

namespace {bool isoCorrect = false;
          
  double ratio_min = .5;
  double ratio_max = 2.5;
           }

void PlotChi2Stat( MnvPlotter &mnvPlotter, MnvH1D* dataHisto, TObjArray mchistos, double dataMCScale, string var, int targetID, int targetZ, bool plotUS, bool plotDS, string playlist, string Name);
void GetMinAndMaxAxisRatio( TH1D* histData, TH1D* histMC, double dataMCScale, double &plotMin, double &plotMax );



//const bool do_fits = true;//Run this code with this boolean as true first and then run it again with boolean as false before plotiing
const bool do_fits = false;//second

TH1D* GetVertErrorBandUniverseHist(MnvH1D* mnvh1d, const std::string& errName, int universe); 

double getPhysChi2( const double * par);
      TH1D* m_histo_trans_data;
      TH1D* m_histo_sig_trans;
      TH1D* m_histo_trans_trans;
      TH1D* m_histo_contin_trans;
      TH1D* m_histo_contin_data;
      TH1D* m_histo_sig_contin;
      TH1D* m_histo_trans_contin;
      TH1D* m_histo_contin_contin;


 
int UnfoldIterations( const std::string& var );


int main(int argc, char * argv[]){
  //ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if(argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n"<<
       "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
      "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" <<
      "\t-Playlist\t \t =\t eg. minervame1A"<< std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    return 0;
  }

  string outdir=argv[1];
  //int targetID = atoi(argv[2]);
  //int targetZ = atoi(argv[3]);
  const string playlist= argv[2];

  TString dbFile;

  dbFile= "/minerva/data/users/zdar/MLPredFiles/ProductionPredictionFiles/me1Dmc/me1Dmc.txt";

 // TFile *fp = new TFile( dbFile,"read" );

  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(playlist);
  NukeCC_Binning  *binsDef = new NukeCC_Binning();
 
  MnvPlotter mnvPlotter(kNukeCCStyle);


   HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(playlist);


   int             id;
   Int_t           run;
   Int_t           subrun;
   Int_t           gate;
   Int_t           slice_number;
   Int_t           ANN_best_segment;
   Float_t         ANN_best_prediction;
   //Int_t           ANN_segments[173];
   //Float_t         ANN_probabilities[173];

   //top 4 segment and probabilities 
   Float_t        ANN_prob_1;
   Int_t           ANN_segment_1;
   Float_t        ANN_prob_2;
   Int_t           ANN_segment_2;
   Float_t        ANN_prob_3;
   Int_t           ANN_segment_3;
   //Float_t        ANN_prob_4;
   //Int_t           ANN_segment_4;

   std::vector<float> ANN_plane_probs;
   std::vector<int> ANN_segments_vec;


   //The input file cern.dat is a copy of the CERN staff data base
   FILE *fp = fopen(dbFile,"r");
   TString filename = Form("%s/MLPred_174Planes_me1Dmc_NX%s.root",outdir.c_str(),playlist.c_str());

   
   cout << "am here 1" << endl;
   TFile *hfile = new TFile(filename, "RECREATE");
   
   TTree *tree = new TTree("MasterAnaDev","MasterAnaDev ML Vertexing Pred");
   
   tree->Branch("run",&run,"run/I");
   tree->Branch("subrun",&subrun,"subrun/I");
   tree->Branch("gate",&gate,"gate/I");
   tree->Branch("slice_number",&slice_number,"slice_number/I");
   //tree->Branch("ANN_probabilities",&ANN_probabilities,"ANN_probabilities/F");
   //tree->Branch("ANN_segments",&ANN_segments,"ANN_segments/I");
   tree->Branch("ANN_best_segment",&ANN_best_segment,"ANN_best_segment/I");
   tree->Branch("ANN_best_prediction",&ANN_best_prediction,"ANN_best_prediction/F");
   tree->Branch("ANN_segment_1",&ANN_segment_1,"ANN_segment_1/I");
   tree->Branch("ANN_prob_1",&ANN_prob_1,"ANN_prob_1/F");
   tree->Branch("ANN_segment_2",&ANN_segment_2,"ANN_segment_2/I");
   tree->Branch("ANN_prob_2",&ANN_prob_2,"ANN_prob_2/F");
   tree->Branch("ANN_segment_3",&ANN_segment_3,"ANN_segment_3/I");
   tree->Branch("ANN_prob_3",&ANN_prob_3,"ANN_prob_3/F");
   
   cout << "am here 2" << endl;
   std::map<float,int> ProbsToSegmentMap;
   float    ANN_probabilities[174];
   
   cout << "am here 3" << endl;
   string last_srunname = "000000_0000";
   char line[5000];
   while (fgets(line,5000,fp)) {
     //sscanf(&line[0],"VALUES(%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f);",
     sscanf(&line[0], "%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", 
     //sscanf(&line[0], "%d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
     &run,&subrun,&gate,&slice_number,&ANN_best_segment,&ANN_probabilities[0],&ANN_probabilities[1],&ANN_probabilities[2],&ANN_probabilities[3],&ANN_probabilities[4],&ANN_probabilities[5],&ANN_probabilities[6],&ANN_probabilities[7],&ANN_probabilities[8],&ANN_probabilities[9],&ANN_probabilities[10],&ANN_probabilities[11],&ANN_probabilities[12],&ANN_probabilities[13],&ANN_probabilities[14],&ANN_probabilities[15],&ANN_probabilities[16],&ANN_probabilities[17],&ANN_probabilities[18],&ANN_probabilities[19],&ANN_probabilities[20],&ANN_probabilities[21],&ANN_probabilities[22],&ANN_probabilities[23],&ANN_probabilities[24],&ANN_probabilities[25],&ANN_probabilities[26],&ANN_probabilities[27],&ANN_probabilities[28],&ANN_probabilities[29],&ANN_probabilities[30],&ANN_probabilities[31],&ANN_probabilities[32],&ANN_probabilities[33],&ANN_probabilities[34],&ANN_probabilities[35],&ANN_probabilities[36],&ANN_probabilities[37],&ANN_probabilities[38],&ANN_probabilities[39],&ANN_probabilities[40],&ANN_probabilities[41],&ANN_probabilities[42],&ANN_probabilities[43],&ANN_probabilities[44],&ANN_probabilities[45],&ANN_probabilities[46],&ANN_probabilities[47],&ANN_probabilities[48],&ANN_probabilities[49],&ANN_probabilities[50],&ANN_probabilities[51],&ANN_probabilities[52],&ANN_probabilities[53],&ANN_probabilities[54],&ANN_probabilities[55],&ANN_probabilities[56],&ANN_probabilities[57],&ANN_probabilities[58],&ANN_probabilities[59],&ANN_probabilities[60],&ANN_probabilities[61],&ANN_probabilities[62],&ANN_probabilities[63],&ANN_probabilities[64],&ANN_probabilities[65],&ANN_probabilities[66],&ANN_probabilities[67],&ANN_probabilities[68],&ANN_probabilities[69],&ANN_probabilities[70],&ANN_probabilities[71],&ANN_probabilities[72],&ANN_probabilities[73],&ANN_probabilities[74],&ANN_probabilities[75],&ANN_probabilities[76],&ANN_probabilities[77],&ANN_probabilities[78],&ANN_probabilities[79],&ANN_probabilities[80],&ANN_probabilities[81],&ANN_probabilities[82],&ANN_probabilities[83],&ANN_probabilities[84],&ANN_probabilities[85],&ANN_probabilities[86],&ANN_probabilities[87],&ANN_probabilities[88],&ANN_probabilities[89],&ANN_probabilities[90],&ANN_probabilities[91],&ANN_probabilities[92],&ANN_probabilities[93],&ANN_probabilities[94],&ANN_probabilities[95],&ANN_probabilities[96],&ANN_probabilities[97],&ANN_probabilities[98],&ANN_probabilities[99],&ANN_probabilities[100],&ANN_probabilities[101],&ANN_probabilities[102],&ANN_probabilities[103],&ANN_probabilities[104],&ANN_probabilities[105],&ANN_probabilities[106],&ANN_probabilities[107],&ANN_probabilities[108],&ANN_probabilities[109],&ANN_probabilities[110],&ANN_probabilities[111],&ANN_probabilities[112],&ANN_probabilities[113],&ANN_probabilities[114],&ANN_probabilities[115],&ANN_probabilities[116],&ANN_probabilities[117],&ANN_probabilities[118],&ANN_probabilities[119],&ANN_probabilities[120],&ANN_probabilities[121],&ANN_probabilities[122],&ANN_probabilities[123],&ANN_probabilities[124],&ANN_probabilities[125],&ANN_probabilities[126],&ANN_probabilities[127],&ANN_probabilities[128],&ANN_probabilities[129],&ANN_probabilities[130],&ANN_probabilities[131],&ANN_probabilities[132],&ANN_probabilities[133],&ANN_probabilities[134],&ANN_probabilities[135],&ANN_probabilities[136],&ANN_probabilities[137],&ANN_probabilities[138],&ANN_probabilities[139],&ANN_probabilities[140],&ANN_probabilities[141],&ANN_probabilities[142],&ANN_probabilities[143],&ANN_probabilities[144],&ANN_probabilities[145],&ANN_probabilities[146],&ANN_probabilities[147],&ANN_probabilities[148],&ANN_probabilities[149],&ANN_probabilities[150],&ANN_probabilities[151],&ANN_probabilities[152],&ANN_probabilities[153],&ANN_probabilities[154],&ANN_probabilities[155],&ANN_probabilities[156],&ANN_probabilities[157],&ANN_probabilities[158],&ANN_probabilities[159],&ANN_probabilities[160],&ANN_probabilities[161],&ANN_probabilities[162],&ANN_probabilities[163],&ANN_probabilities[164],&ANN_probabilities[165],&ANN_probabilities[166],&ANN_probabilities[167],&ANN_probabilities[168],&ANN_probabilities[169],&ANN_probabilities[170],&ANN_probabilities[171],&ANN_probabilities[172],&ANN_probabilities[173]);
     
     ProbsToSegmentMap.clear();
     ANN_best_prediction = ANN_probabilities[ANN_best_segment];
     
     int ANN_segments[174];
     for(int m=0; m<174; m++ ){
       ANN_segments[m] = m;
       ProbsToSegmentMap[ANN_probabilities[m]] = m;
     }
     
     /* for checking hdf5 entries
	string srunname = Form("%d_%d",run,subrun);
	if( srunname != last_srunname ){
	cout << srunname << endl;
	last_srunname = srunname;
	}
	else continue;
     */
     
     int index = 0; 
     ANN_segments_vec.clear();
     ANN_plane_probs.clear();
     
     for( std::map<float,int>::reverse_iterator rit = ProbsToSegmentMap.rbegin(); rit != ProbsToSegmentMap.rend(); ++rit ){
       ANN_plane_probs.push_back( rit->first );
       ANN_segments_vec.push_back( rit->second );
       index++;
       if( index > 3 ) break; 
     }
     
     ANN_segment_1 = ANN_segments_vec[0];
     ANN_prob_1 = ANN_plane_probs[0];
     ANN_segment_2 = ANN_segments_vec[1];
     ANN_prob_2 = ANN_plane_probs[1];
     ANN_segment_3 = ANN_segments_vec[2];
     ANN_prob_3 = ANN_plane_probs[2];
     
     
       cout << "run, subrun, gate, slice = " << run << ", " << subrun << ", " << gate << ", " << slice_number;
 //      cout << "     best segment, probs = " << ANN_best_segment << ", " << ANN_best_prediction << endl;
 //      cout << "     segment 1, prob 1 = " << ANN_segment_1 << ", " << ANN_prob_1 << endl;
 //      cout << "     vector segments 0, prob segments 0 = " << ANN_segments_vec[0] << ", " << ANN_plane_probs[0] << endl;
       for(int m=0; m<174; m++ )
         cout << "segment, probs = " << m << ", " << ANN_probabilities[m] << endl;  
     
     tree->Fill();
   }
   tree->Write();
   
   fclose(fp);
   delete tree;
   delete hfile;
   return 0;

}


