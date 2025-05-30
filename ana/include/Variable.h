#ifndef VARIABLE_H
#define VARIABLE_H

#include <iterator>
#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "PlotUtils/HistFolio.h"
#include "PlotUtils/HistWrapper.h"

#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MnvH2D.h"
#ifndef __CINT__  // CINT doesn't know about std::function
#include "PlotUtils/VariableBase.h"
#include "PlotUtils/Variable2DBase.h"
#include "MinervaUnfold/MnvResponse.h"
#include <bitset>
#include "PlotUtils/AnaBinning.h"

#endif  // __CINT__

namespace VarLoop {

class Variable : public PlotUtils::VariableBase<NUKECC_ANA::CVUniverse> {
 private:
  typedef PlotUtils::HistWrapper<NUKECC_ANA::CVUniverse> HW;
  typedef PlotUtils::MnvH1D MH1D;

 public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  Variable(ARGS... args) : PlotUtils::VariableBase<NUKECC_ANA::CVUniverse>(args...) {}

  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // HISTWRAPPER
  // selected mc reco histwrapper
  HW m_selected_mc_reco, m_selected_data_reco, m_selected_truth_reco, m_selected_data_reco_sb, m_selected_mc_true, m_selected_mc_reco_USCH, m_selected_mc_reco_DSCH, m_selected_mc_reco_other, m_selected_trans_SB, m_selected_contin_SB, m_selected_mc_reco_true_target;

  typedef PlotUtils::Hist2DWrapper<NUKECC_ANA::CVUniverse> HW2D;
  HW2D mresp1D;
  HW2D m_selected_Migration;
  typedef PlotUtils::MnvH2D MH2D;
  std::map<std::string,MinervaUnfold::MnvResponse*> Response1D;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr2;

  MnvH2D *migrationH2D = NULL;
  MnvH1D *h_reco1D = NULL;
  MnvH1D *h_truth1D = NULL;
 // HISTFOLIO
  // selected mc reco - signal background histfolio
  PlotUtils::HistFolio<MH1D> m_selected_mc_sb;
  // PlotUtils::MH1D* m_selected_data_sb;
  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHistograms(T univs) {
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    const bool clear_bands = true;  // we want empty histograms
    std::cout<<"Before"<<std::endl;
    // HISTWRAPPER
    // For background from runEventloop script selected
    //MH1D* dummy_selected_mc_reco_USCH = new MH1D(Form("h_USCH_%s", name), name,GetNBins(), bins.data());
    //m_selected_mc_reco_USCH = HW(dummy_selected_mc_reco_USCH, univs, clear_bands);

    // For background from runEventloop script selected
    //MH1D* dummy_selected_mc_reco_DSCH = new MH1D(Form("h_DSCH_%s", name), name,GetNBins(), bins.data());
    //m_selected_mc_reco_DSCH = HW(dummy_selected_mc_reco_DSCH, univs, clear_bands);

    // For background from runEventloop script selected
    //MH1D* dummy_selected_mc_reco_other = new MH1D(Form("h_other_%s", name), name,GetNBins(), bins.data());
    //m_selected_mc_reco_other = HW(dummy_selected_mc_reco_other, univs, clear_bands);
    std::cout<<"Before 1"<<std::endl;
    // For background from runEventloop script selected
    //MH1D* dummy_selected_hist_true_trans_SB = new MH1D(Form("true_trans_in_target_%s", name), name,GetNBins(), bins.data());
    //m_selected_trans_SB = HW(dummy_selected_hist_true_trans_SB, univs, clear_bands);

    // For background from runEventloop script selected
    //MH1D* dummy_selected_hist_true_contin_SB = new MH1D(Form("true_contin_in_target_%s", name), name,GetNBins(), bins.data());
    //m_selected_contin_SB = HW(dummy_selected_hist_true_contin_SB, univs, clear_bands);

    // For background from runEventloop script selected
    MH1D* dummy_selected_mc_reco_true_target = new MH1D(Form("true_target_%s", name), name,GetNBins(), bins.data());
    m_selected_mc_reco_true_target = HW(dummy_selected_mc_reco_true_target, univs, clear_bands);
    std::cout<<"Before 2"<<std::endl;
    // selected mc reco histwrapper
    MH1D* dummy_selected_mc_reco = new MH1D(Form("h_mc_%s", name), name,GetNBins(), bins.data());
    m_selected_mc_reco = HW(dummy_selected_mc_reco, univs, clear_bands);
    std::cout<<"Before 3"<<std::endl;
    //for 1D analysis
    MH1D* dummy_selected_mc_true = new MH1D(Form("h_mc_true_%s", name), name,GetNBins(), bins.data());
    m_selected_mc_true = HW(dummy_selected_mc_true, univs, clear_bands);
    std::cout<<"Before 4"<<std::endl;
    //////////For Efficiency denominator
    MH1D* dummy_selected_truth_reco = new MH1D(Form("h_truth_%s", name), name, GetNBins(), bins.data());
    m_selected_truth_reco = HW(dummy_selected_truth_reco, univs, clear_bands);
    std::cout<<"Before 5"<<std::endl;

    //For Data
    MH1D* dummy_selected_data_reco = new MH1D(Form("h_data_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco = HW(dummy_selected_data_reco, univs, clear_bands);
    std::cout<<"Before 6"<<std::endl;

    //For Data in sidebands
    //MH1D* selected_data_reco_sb = new MH1D(Form("h_data_sb_%s", name), name, GetNBins(), bins.data());
    //m_selected_data_reco_sb = HW(selected_data_reco_sb, univs, clear_bands);

    //HISTFOLIO
    //selected mc reco - signal background histfolio
    m_selected_mc_sb = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("selected_mc_sb_%s", name), name, GetNBins(), bins.data());

 //PlotUtils::MnvH1D* data = new PlotUtils::MnvH1D(
 //   "dummy", "dummy", plotting::nbins, plotting::xmin, plotting::xmax);
 //  m_selected_data_sb = PlotUtils::HistFolio<PlotUtils::MnvH1D>(
 //    Form("selected_data_sb_%s", name), name, GetNBins(), bins.data());

    m_selected_mc_sb.AddComponentHist("DIS");
    m_selected_mc_sb.AddComponentHist("MC");
   // m_selected_data_sb.AddComponentHist("Data");
   
   MH2D* dummy_selected_Migration = new MH2D(Form("selected_Migration_%s", name), name, GetNBins(),GetBinVec().data(), GetNBins(), GetBinVec().data());
    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);

delete dummy_selected_mc_reco;
delete dummy_selected_truth_reco;
delete dummy_selected_mc_true;
delete dummy_selected_data_reco;
delete dummy_selected_Migration;
std::cout<<"Hola Oscar"<<std::endl;
  }

//=======================================================================================
  // When plotting, don't make Cuts on this Variable's axes
  //=======================================================================================
  //TODO: Move into PlotUtils?
  #define cut_map_t std::bitset<64>
  std::bitset<64> ignoreTheseCuts;
  //Populate ignoreTheseCuts.
  template <class CUTS_T> //Kill two bird templates with one stone template
  cut_map_t MatchCutsToVars(const CUTS_T& cuts)
  {
     const char* name = GetName().c_str();
    //Stupid compiler!  Let me use generic lambdas!
    using cut_t = decltype(*cuts.begin());
    ignoreTheseCuts.reset(); //Make sure I start off ignoring NO Cuts.
                             //Only cuts with 1 will be ignored.
    const auto found = std::find_if(cuts.begin(), cuts.end(),
                                    [this](const cut_t cut)
                                    {
                                      return cut->getName() == this->GetName();
                                    });

    if(found != cuts.end()) //If we found a match
    {
      ignoreTheseCuts.flip(std::distance(cuts.begin(), found));
    }
    else std::cerr << "Failed to find ANY Cut that matches " << GetName() << ".  You should double-check your Cut names.\n";
    return ignoreTheseCuts; //Just here in case it helps with debugging
  }
  //Use ignoreTheseCuts
  inline cut_map_t IgnoreMyVars(const cut_map_t allCuts) const
  {
    return allCuts | ignoreTheseCuts;
  }




//void SetupResponse1D(std::map<const std::string, const int> systematics){
void SetupResponse1D(std::map<const std::string, int> systematics){
//void SetupResponse(T univs){
	   std::vector<double> bins = GetBinVec();
	   const char* name = GetName().c_str();

	   axis_binning bin_x;
	   bin_x.uniform=false;

     	   vector<double> vx;

	   for(int i=0; i<=GetNBins(); i++){vx.push_back(GetBinVec().data()[i]);}

	   bin_x.bin_edges = vx;
	   bin_x.nbins	   = GetNBins();
	   bin_x.min 	   = GetBinVec().data()[0];
	   bin_x.max       = GetBinVec().data()[GetNBins()];

           cout<<"bins min = "<<bin_x.min<<"bins max = "<<bin_x.max<<"Number of bins = "<<bin_x.nbins<<endl;

	   //Response1D.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("selected_mc_response1d_%s", name), name, bin_x, bin_x, systematics)));
	   Response1D.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("response1d_%s", name), name, bin_x, bin_x, systematics)));
}


//===================================================================================
void FillResponse1D(double x_reco, double x_true, const std::string name, double w, int unv){

	//std::cout<<name<<std::endl;
 	for(mnv_itr = Response1D.begin(); mnv_itr != Response1D.end(); ++mnv_itr){
		(mnv_itr->second)->Fill(x_reco,x_true,name,unv, w);
	}
}

//=====================================

template <typename T>
void getResponseObjects1D(T univs)
{
//  bool status = false;
  for(mnv_itr2 = Response1D.begin(); mnv_itr2 != Response1D.end(); ++mnv_itr2){
                (mnv_itr2->second)->GetMigrationObjects( migrationH2D, h_reco1D, h_truth1D );
        }
  const bool clear_bands = true;
  mresp1D = HW2D(migrationH2D, univs, clear_bands);
}

  //=======================================================================================
  // WRITE ALL HISTOGRAMS FOR EVENTLOOP
  //=======================================================================================
 void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC) { m_selected_mc_reco.hist->Write();
               //m_selected_mc_reco_USCH.hist->Write();
               //m_selected_mc_reco_DSCH.hist->Write();
               //m_selected_mc_reco_other.hist->Write();
               //m_selected_trans_SB.hist->Write();
               //m_selected_contin_SB.hist->Write();
               m_selected_mc_reco_true_target.hist->Write();
}
    else {m_selected_data_reco.hist->Write();
}
    // selected mc  histfolio fir Hist Stacking
   //if(isMC) m_selected_mc_sb.WriteToFile(f);
   //else m_selected_data_reco_sb.hist->Write();
  }
  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
 void WriteAllHistogramsToFileEff(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC) { m_selected_mc_reco.hist->Write();
    }
    else { m_selected_truth_reco.hist->Write();
    }
}
  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
 void WriteAllHistogramsToFileMig(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC) { m_selected_mc_reco.hist->Write();
               m_selected_mc_true.hist->Write();
               mresp1D.hist->Write();
}
    else {m_selected_data_reco.hist->Write();
}
    // selected mc  histfolio fir Hist Stacking
   //if(isMC) m_selected_mc_sb.WriteToFile(f);
   //else m_selected_data_reco_sb.hist->Write();
  }
};



}  // namespace VarLoop



namespace Var2DLoop {
class Variable2D : public PlotUtils::Variable2DBase<NUKECC_ANA::CVUniverse> {
 private:
  //=======================================================================================
  // TYPEDEFS CONVENIENCE
  //=======================================================================================
  typedef PlotUtils::Hist2DWrapper<NUKECC_ANA::CVUniverse> HW2D;
  //typedef Histogram<NUKECC_ANA::CVUniverse> HW2D;
  typedef PlotUtils::MnvH2D MH2D;

 public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  Variable2D(ARGS... args) : Variable2DBase<NUKECC_ANA::CVUniverse>(args...) {}

  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // HISTWRAPPER
  HW2D m_selected_mc_reco, m_selected_mc_reco_mig, m_selected_mc_true_mig,  m_selected_data_reco,m_selected_truth_reco,m_selected_Migration, mresp;
  HW2D m_selected_mc_reco_true_target, m_selected_mc_reco_USCH, m_selected_mc_reco_DSCH, m_selected_mc_reco_other, m_selected_trans_SB, m_selected_contin_SB;
  HW2D m_selected_mc_reco_Lead,m_selected_mc_reco_Carbon,m_selected_mc_reco_Other;

  HW2D m_selected_mc_reco_QE, m_selected_mc_reco_RES, m_selected_mc_reco_DIS, m_selected_mc_reco_2p2h, m_selected_mc_reco_OtherIT;
  HW2D m_selected_mc_truth_QE, m_selected_mc_truth_RES, m_selected_mc_truth_DIS, m_selected_mc_truth_2p2h, m_selected_mc_truth_OtherIT;

  HW2D m_selected_mc_reco_NoScattering, m_selected_mc_reco_ChargeExt, m_selected_mc_reco_Elasticity, m_selected_mc_reco_Absorption;
  HW2D m_selected_mc_reco_PionProd, m_selected_mc_reco_MultNuc, m_selected_mc_reco_KnockOut;
  HW2D m_selected_mc_reco_fate_nn_50, m_selected_mc_reco_fate_pn_51, m_selected_mc_reco_fate_pp_52;

  HW2D m_selected_mc_truth_reco_NoScattering, m_selected_mc_truth_reco_ChargeExt, m_selected_mc_truth_reco_Elasticity, m_selected_mc_truth_reco_Absorption;
  HW2D m_selected_mc_truth_reco_PionProd,     m_selected_mc_truth_reco_MultNuc,   m_selected_mc_truth_reco_KnockOut;
  HW2D m_selected_mc_truth_reco_fate_nn_50, m_selected_mc_truth_reco_fate_pn_51, m_selected_mc_truth_reco_fate_pp_52;
  //MinervaUnfold::MnvResponse* Response;
  std::map<std::string,MinervaUnfold::MnvResponse*> Response;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr2;

 
  MnvH2D *migrationH = NULL;
  MnvH2D *h_reco = NULL;
  MnvH2D *h_truth = NULL;
  //// HISTFOLIO
  // PlotUtils::HistFolio<MH2D> m_selected_mc_sb;

  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHistograms(T univs) {
    const bool clear_bands = true;  // we want empty histograms
    const char* name = GetName().c_str();

   //SideBand histograms for 2D analysis
    MH2D* dummy_selected_mc_reco_true_target =
        new MH2D(Form("h_mc_true_target_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_true_target = HW2D(dummy_selected_mc_reco_true_target, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_USCH =
        new MH2D(Form("h_mc_USCH_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_USCH = HW2D(dummy_selected_mc_reco_USCH, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_DSCH =
        new MH2D(Form("h_mc_DSCH_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_DSCH = HW2D(dummy_selected_mc_reco_DSCH, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_other =
        new MH2D(Form("h_mc_other_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_other = HW2D(dummy_selected_mc_reco_other, univs, clear_bands);

    MH2D* dummy_selected_trans_SB =
        new MH2D(Form("h_mc_trans_SB_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_trans_SB = HW2D(dummy_selected_trans_SB, univs, clear_bands);

    MH2D* dummy_selected_contin_SB =
        new MH2D(Form("h_mc_contin_SB_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_contin_SB = HW2D(dummy_selected_contin_SB, univs, clear_bands);

    // HISTWRAPPER
    // selected mc reco histwrapper
    MH2D* dummy_selected_mc_reco =
        new MH2D(Form("h_mc_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco = HW2D(dummy_selected_mc_reco, univs, clear_bands);

    //For reco MC in reco space
    MH2D* dummy_selected_mc_reco_mig =
        new MH2D(Form("h_mc_reco_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_mig = HW2D(dummy_selected_mc_reco_mig, univs, clear_bands);
    //For reco MC in true space
    MH2D* dummy_selected_mc_true_mig =
        new MH2D(Form("h_mc_true_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_true_mig = HW2D(dummy_selected_mc_true_mig, univs, clear_bands);

    MH2D* dummy_selected_truth_reco =
        new MH2D(Form("h_truth_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_truth_reco = HW2D(dummy_selected_truth_reco, univs, clear_bands);

    //For reco MC in reco for Lead events=======================================================================================
    MH2D* dummy_selected_mc_reco_Lead =
        new MH2D(Form("h_mc_reco_Lead_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_Lead = HW2D(dummy_selected_mc_reco_Lead, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_Carbon =
        new MH2D(Form("h_mc_reco_Carbon_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_Carbon = HW2D(dummy_selected_mc_reco_Carbon, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_Other =
        new MH2D(Form("h_mc_reco_Other_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_Other = HW2D(dummy_selected_mc_reco_Other, univs, clear_bands);

    //Channel Interaction=======================================================================================================
    MH2D* dummy_m_selected_mc_reco_QE =
        new MH2D(Form("m_selected_mc_reco_QE_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_QE = HW2D(dummy_m_selected_mc_reco_QE, univs, clear_bands);

    MH2D* dummy_m_selected_mc_reco_RES =
        new MH2D(Form("m_selected_mc_reco_RES_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_RES = HW2D(dummy_m_selected_mc_reco_RES, univs, clear_bands);

    MH2D* dummy_m_selected_mc_reco_DIS =
        new MH2D(Form("m_selected_mc_reco_DIS_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_DIS = HW2D(dummy_m_selected_mc_reco_DIS, univs, clear_bands);

    MH2D* dummy_m_selected_mc_reco_2p2h =
        new MH2D(Form("m_selected_mc_reco_2p2h_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_2p2h = HW2D(dummy_m_selected_mc_reco_2p2h, univs, clear_bands);

    MH2D* dummy_m_selected_mc_reco_OtherIT =
        new MH2D(Form("m_selected_mc_reco_OtherIT_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_OtherIT = HW2D(dummy_m_selected_mc_reco_OtherIT, univs, clear_bands);

    MH2D* dummy_m_selected_mc_truth_QE =
        new MH2D(Form("m_selected_mc_truth_QE_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_QE = HW2D(dummy_m_selected_mc_truth_QE, univs, clear_bands);

    MH2D* dummy_m_selected_mc_truth_RES =
        new MH2D(Form("m_selected_mc_truth_RES_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_RES = HW2D(dummy_m_selected_mc_truth_RES, univs, clear_bands);

    MH2D* dummy_m_selected_mc_truth_DIS =
        new MH2D(Form("m_selected_mc_truth_DIS_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_DIS = HW2D(dummy_m_selected_mc_truth_DIS, univs, clear_bands);

    MH2D* dummy_m_selected_mc_truth_2p2h =
        new MH2D(Form("m_selected_mc_truth_2p2h_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_2p2h = HW2D(dummy_m_selected_mc_truth_2p2h, univs, clear_bands);

    MH2D* dummy_m_selected_mc_truth_OtherIT =
        new MH2D(Form("m_selected_mc_truth_OtherIT_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_OtherIT = HW2D(dummy_m_selected_mc_truth_OtherIT, univs, clear_bands);

    //For reco MC in reco for Fate breakdown=======================================================================================
    MH2D* dummy_selected_mc_reco_NoScattering =
        new MH2D(Form("h_mc_reco_NoScattering_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_NoScattering = HW2D(dummy_selected_mc_reco_NoScattering, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_ChargeExt =
        new MH2D(Form("h_mc_reco_ChargeExt_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_ChargeExt = HW2D(dummy_selected_mc_reco_ChargeExt, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_Elasticity =
        new MH2D(Form("h_mc_reco_Elasticity_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_Elasticity = HW2D(dummy_selected_mc_reco_Elasticity, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_Absorption =
        new MH2D(Form("h_mc_reco_Absorption_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_Absorption = HW2D(dummy_selected_mc_reco_Absorption, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_PionProd =
        new MH2D(Form("h_mc_reco_PionProd_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_PionProd = HW2D(dummy_selected_mc_reco_PionProd, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_MultNuc =
        new MH2D(Form("h_mc_reco_MultNuc_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_MultNuc = HW2D(dummy_selected_mc_reco_MultNuc, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_KnockOut =
        new MH2D(Form("h_mc_reco_KnockOut_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_KnockOut = HW2D(dummy_selected_mc_reco_KnockOut, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_fate_nn_50 =
        new MH2D(Form("h_mc_reco_fate_nn_50_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_fate_nn_50 = HW2D(dummy_selected_mc_reco_fate_nn_50, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_fate_pn_51 =
        new MH2D(Form("h_mc_reco_fate_pn_51_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_fate_pn_51 = HW2D(dummy_selected_mc_reco_fate_pn_51, univs, clear_bands);

    MH2D* dummy_selected_mc_reco_fate_pp_52 =
        new MH2D(Form("h_mc_reco_fate_pp_52_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco_fate_pp_52 = HW2D(dummy_selected_mc_reco_fate_pp_52, univs, clear_bands);

    //For reco MC in reco for Fate breakdown======================================================================================= 
    MH2D* dummy_selected_mc_truth_reco_NoScattering =
        new MH2D(Form("h_mc_truth_reco_NoScattering_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_reco_NoScattering = HW2D(dummy_selected_mc_truth_reco_NoScattering, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_reco_ChargeExt =
        new MH2D(Form("h_mc_truth_reco_ChargeExt_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_reco_ChargeExt = HW2D(dummy_selected_mc_truth_reco_ChargeExt, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_reco_Elasticity =
        new MH2D(Form("h_mc_truth_reco_Elasticity_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_reco_Elasticity = HW2D(dummy_selected_mc_truth_reco_Elasticity, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_reco_Absorption =
        new MH2D(Form("h_mc_truth_reco_Absorption_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_reco_Absorption = HW2D(dummy_selected_mc_truth_reco_Absorption, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_reco_PionProd =
        new MH2D(Form("h_mc_truth_reco_PionProd_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_reco_PionProd = HW2D(dummy_selected_mc_truth_reco_PionProd, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_reco_MultNuc =
        new MH2D(Form("h_mc_truth_reco_MultNuc_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_reco_MultNuc = HW2D(dummy_selected_mc_truth_reco_MultNuc, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_reco_KnockOut =
        new MH2D(Form("h_mc_truth_reco_KnockOut_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_reco_KnockOut = HW2D(dummy_selected_mc_truth_reco_KnockOut, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_reco_fate_nn_50 =
        new MH2D(Form("h_mc_truth_reco_fate_nn_50_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_reco_fate_nn_50 = HW2D(dummy_selected_mc_truth_reco_fate_nn_50, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_reco_fate_pn_51 =
        new MH2D(Form("h_mc_truth_reco_fate_pn_51_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_reco_fate_pn_51 = HW2D(dummy_selected_mc_truth_reco_fate_pn_51, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_reco_fate_pp_52 =
        new MH2D(Form("h_mc_truth_reco_fate_pp_52_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_reco_fate_pp_52 = HW2D(dummy_selected_mc_truth_reco_fate_pp_52, univs, clear_bands);

//===============================================================================================================================
   // MH2D* dummy_selected_truth_reco =
   //     new MH2D(Form("selected_truth2d_reco_%s", name), name, setBinLogX(), GetBinVecX().data(), setBinLogY(), GetBinVecY().data());
   //              GetBinVecX().data(), setBinLogY(), GetBinVecY().data());
   // m_selected_truth_reco = HW2D(dummy_selected_truth_reco, univs, clear_bands);

    MH2D* dummy_selected_data_reco =
        new MH2D(Form("h_data_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_data_reco = HW2D(dummy_selected_data_reco, univs, clear_bands);


    //MH2D* dummy_selected_Migration =
   //     new MH2D(Form("selected_Migration_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
   // m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);


    /////For 2D analysis
    MH2D* dummy_selected_Migration = new MH2D(Form("selected_Migration_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);///For 2D Migration


    ///For 1D analysis
    /*MH2D* dummy_selected_Migration =
        new MH2D(Form("selected_Migration_%s", name), name, GetNBinsX(),GetBinVecX().data(), GetNBinsX(), GetBinVecX().data());
    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);*/

    //delete Response;
    delete dummy_selected_mc_reco;
    delete dummy_selected_truth_reco;
    delete dummy_selected_Migration;
    delete dummy_selected_data_reco;
    delete dummy_selected_mc_reco_true_target;
    delete dummy_selected_mc_reco_USCH;
    delete dummy_selected_mc_reco_DSCH;
    delete dummy_selected_mc_reco_other;
    delete dummy_selected_trans_SB;
    delete dummy_selected_contin_SB;
    delete dummy_selected_mc_reco_mig;
    delete dummy_selected_mc_true_mig;

    delete dummy_m_selected_mc_reco_QE;
    delete dummy_m_selected_mc_reco_RES;
    delete dummy_m_selected_mc_reco_DIS;
    delete dummy_m_selected_mc_reco_2p2h;
    delete dummy_m_selected_mc_reco_OtherIT;

    delete dummy_m_selected_mc_truth_QE;
    delete dummy_m_selected_mc_truth_RES;
    delete dummy_m_selected_mc_truth_DIS;
    delete dummy_m_selected_mc_truth_2p2h;
    delete dummy_m_selected_mc_truth_OtherIT;

    delete dummy_selected_mc_reco_NoScattering;
    delete dummy_selected_mc_reco_ChargeExt;
    delete dummy_selected_mc_reco_Elasticity;
    delete dummy_selected_mc_reco_Absorption;
    delete dummy_selected_mc_reco_PionProd;
    delete dummy_selected_mc_reco_MultNuc;
    delete dummy_selected_mc_reco_KnockOut;
    delete dummy_selected_mc_reco_fate_nn_50;
    delete dummy_selected_mc_reco_fate_pn_51;
    delete dummy_selected_mc_reco_fate_pp_52;

    delete dummy_selected_mc_truth_reco_NoScattering;
    delete dummy_selected_mc_truth_reco_ChargeExt;
    delete dummy_selected_mc_truth_reco_Elasticity;
    delete dummy_selected_mc_truth_reco_Absorption;
    delete dummy_selected_mc_truth_reco_PionProd;
    delete dummy_selected_mc_truth_reco_MultNuc;
    delete dummy_selected_mc_truth_reco_KnockOut;
    delete dummy_selected_mc_truth_reco_fate_nn_50;
    delete dummy_selected_mc_truth_reco_fate_pn_51;
    delete dummy_selected_mc_truth_reco_fate_pp_52;

  }
//=====



//void SetupResponse(std::map<const std::string, const int> systematics){
//After the change the below one was included
void SetupResponse(std::map<const std::string, int> systematics){
//void SetupResponse(T univs){

	   const char* name = GetName().c_str();
	   axis_binning bin_x, bin_y;
	   bin_x.uniform=false;

	   vector<double> vx;

	   for(int i=0; i<=GetNBinsX(); i++){vx.push_back(GetBinVecX().data()[i]);}

	   vector<double> vy;
	   for(int j=0; j<=GetNBinsY(); j++){vy.push_back(GetBinVecY().data()[j]);}
	   bin_x.bin_edges = vx;
	   bin_x.nbins	    = GetNBinsX();
	   bin_x.min 	    = GetBinVecX().data()[0];
	   bin_x.max       = GetBinVecX().data()[GetNBinsX()];
	   bin_y.bin_edges = vy;
	   bin_y.nbins	    = GetNBinsY();
	   bin_y.min 	    = GetBinVecY().data()[0];
	   bin_y.max       = GetBinVecY().data()[GetNBinsY()];

std::cout << "vx.size(): " << vx.size() << ", bin_x.nbins: " << bin_x.nbins << std::endl;
for (size_t i = 1; i < vx.size(); ++i) {
    if (vx[i] <= vx[i - 1])
        std::cerr << "ERROR: Bin edge not strictly increasing at " << i << std::endl;
}



	  std::cout << "bin_x edges:\n";
for (auto& b : bin_x.bin_edges) std::cout << b << " ";
std::cout << "\n";

std::cout << "bin_y edges:\n";
for (auto& b : bin_y.bin_edges) std::cout << b << " ";
std::cout << "\n";


	   //Response.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("selected_mc_response2d_%s", name), name, bin_x, bin_y, bin_x, bin_y, systematics)));
	   Response.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("selected_mc_response2d_%s", name), name, bin_x, bin_y, bin_x, bin_y, systematics)));
	   //Response.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("response2d_%s", name), name, bin_x, bin_y, bin_x, bin_y, systematics)));
}

//===================================================================================
//
//===================================================================================
void FillResponse(double x_reco, double y_reco, double x_true, double y_true, const std::string name, double w, int unv){
 	for(mnv_itr = Response.begin(); mnv_itr != Response.end(); ++mnv_itr){
		(mnv_itr->second)->Fill(x_reco,y_reco,x_true,y_true,name,unv, w);
	}


}
//===================================================================================
//
//===================================================================================
template <typename T>
void getResponseObjects(T univs)
{
 // bool status = false;
  for(mnv_itr2 = Response.begin(); mnv_itr2 != Response.end(); ++mnv_itr2){
                (mnv_itr2->second)->GetMigrationObjects( migrationH, h_reco, h_truth );;
        }
  const bool clear_bands = true;
  mresp = HW2D(migrationH, univs, clear_bands);
}

//=======================================================================================
// WRITE ALL HISTOGRAMS
//=======================================================================================
void WriteAllHistogramsToFile(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ m_selected_mc_reco.hist->Write();
                 m_selected_mc_reco_Lead.hist->Write();
                 m_selected_mc_reco_Carbon.hist->Write();
                 m_selected_mc_reco_Other.hist->Write();
                 m_selected_mc_reco_true_target.hist->Write();
                 m_selected_mc_reco_NoScattering.hist->Write();
                 m_selected_mc_reco_ChargeExt.hist->Write();
                 m_selected_mc_reco_Elasticity.hist->Write();
                 m_selected_mc_reco_Absorption.hist->Write();
                 m_selected_mc_reco_PionProd.hist->Write();
                 m_selected_mc_reco_MultNuc.hist->Write();
                 m_selected_mc_reco_KnockOut.hist->Write();
                 m_selected_mc_reco_fate_nn_50.hist->Write();
                 m_selected_mc_reco_fate_pn_51.hist->Write();
                 m_selected_mc_reco_fate_pp_52.hist->Write();

                 m_selected_mc_truth_reco_NoScattering.hist->Write();
                 m_selected_mc_truth_reco_ChargeExt.hist->Write();
                 m_selected_mc_truth_reco_Elasticity.hist->Write();
                 m_selected_mc_truth_reco_Absorption.hist->Write();
                 m_selected_mc_truth_reco_PionProd.hist->Write();
                 m_selected_mc_truth_reco_MultNuc.hist->Write();
                 m_selected_mc_truth_reco_KnockOut.hist->Write();
                 m_selected_mc_truth_reco_fate_nn_50.hist->Write();
                 m_selected_mc_truth_reco_fate_pn_51.hist->Write();
                 m_selected_mc_truth_reco_fate_pp_52.hist->Write();

                 m_selected_mc_reco_QE.hist->Write();
                 m_selected_mc_reco_RES.hist->Write();
                 m_selected_mc_reco_DIS.hist->Write();
                 m_selected_mc_reco_2p2h.hist->Write();
                 m_selected_mc_reco_OtherIT.hist->Write();

                 m_selected_mc_truth_QE.hist->Write();
                 m_selected_mc_truth_RES.hist->Write();
                 m_selected_mc_truth_DIS.hist->Write();
                 m_selected_mc_truth_2p2h.hist->Write();
                 m_selected_mc_truth_OtherIT.hist->Write();
                 //m_selected_mc_reco_USCH.hist->Write();
                 //m_selected_mc_reco_DSCH.hist->Write();
                 //m_selected_mc_reco_other.hist->Write();
                 //m_selected_trans_SB.hist->Write();
                 //m_selected_contin_SB.hist->Write();
   }//  Response->GetMigrationMatrix()->Write();}
    else { m_selected_data_reco.hist->Write();
     }
    // selected mc reco
  }
//=======================================================================================
// WRITE ALL HISTOGRAMS
//=======================================================================================
void WriteAllHistogramsToFileEff(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ m_selected_mc_reco.hist->Write();
   }
    else { m_selected_truth_reco.hist->Write();
     }
    // selected mc reco
  }
//=======================================================================================
// WRITE ALL HISTOGRAMS
//=======================================================================================
void WriteAllHistogramsToFileMig(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ m_selected_mc_reco.hist->Write();
                 m_selected_Migration.hist->Write();
                 m_selected_mc_reco_mig.hist->Write();
                 m_selected_mc_true_mig.hist->Write();
      //for (auto responses:Response)responses.second->GetMigrationMatrix()->Write();
        mresp.hist->Write();
   }//  Response->GetMigrationMatrix()->Write();}
    else {m_selected_data_reco.hist->Write();
     }
  }
};
}  // namespace Var2DLoop


#endif  // VARIABLE_H
