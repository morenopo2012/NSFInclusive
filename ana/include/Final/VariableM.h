#ifndef VARIABLE_H
#define VARIABLE_H

#include "include/CVUniverse.h"
#include "PlotUtils/HistFolio.h"
#include "PlotUtils/HistWrapper.h"

#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MnvH2D.h"
#ifndef __CINT__  // CINT doesn't know about std::function
#include "PlotUtils/VariableBase.h"

#include "PlotUtils/Variable2DBase.h"
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
  HW m_selected_mc_reco,m_selected_data_reco, m_selected_truth_reco, m_selected_data_reco_sb;

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

    // HISTWRAPPER
    // selected mc reco histwrapper
    MH1D* dummy_selected_mc_reco = new MH1D(Form("selected_mc_reco_%s", name), name,
                                      GetNBins(), bins.data());
    m_selected_mc_reco = HW(dummy_selected_mc_reco, univs, clear_bands);

    MH1D* dummy_selected_truth_reco = new MH1D(Form("selected_truth_reco_%s", name), name, GetNBins(), bins.data());
    m_selected_truth_reco = HW(dummy_selected_truth_reco, univs, clear_bands);
    
    
    MH1D* dummy_selected_data_reco = new MH1D(Form("selected_data_reco_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco = HW(dummy_selected_data_reco, univs, clear_bands);


    MH1D* selected_data_reco_sb = new MH1D(Form("selected_data_reco_sb_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco_sb = HW(selected_data_reco_sb, univs, clear_bands);
  
  // HISTFOLIO
    // selected mc reco - signal background histfolio
    
   m_selected_mc_sb = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("selected_mc_sb_%s", name), name, GetNBins(), bins.data());
    
 // PlotUtils::MnvH1D* data = new PlotUtils::MnvH1D(
   //   "dummy", "dummy", plotting::nbins, plotting::xmin, plotting::xmax);
  //  m_selected_data_sb = PlotUtils::HistFolio<PlotUtils::MnvH1D>(
    //    Form("selected_data_sb_%s", name), name, GetNBins(), bins.data());
   
    m_selected_mc_sb.AddComponentHist("DIS");
    m_selected_mc_sb.AddComponentHist("MC");
   // m_selected_data_sb.AddComponentHist("Data");
delete dummy_selected_mc_reco;
delete dummy_selected_truth_reco;
delete dummy_selected_data_reco;
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
 void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC)  m_selected_mc_reco.hist->Write();
    else m_selected_truth_reco.hist->Write();

    // selected mc  histfolio fir Hist Stacking
   if(isMC) m_selected_mc_sb.WriteToFile(f);
   else m_selected_data_reco_sb.hist->Write();
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
  HW2D m_selected_mc_reco,m_selected_data_reco, m_selected_truth_reco;

  //// HISTFOLIO
  // PlotUtils::HistFolio<MH2D> m_selected_mc_sb;

  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHistograms(T univs) {
    const bool clear_bands = true;  // we want empty histograms
    const char* name = GetName().c_str();

    // HISTWRAPPER
    // selected mc reco histwrapper
    MH2D* dummy_selected_mc_reco =
        new MH2D(Form("selected_mc_reco2d_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco = HW2D(dummy_selected_mc_reco, univs, clear_bands);
    
    
    MH2D* dummy_selected_truth_reco =
        new MH2D(Form("selected_truth2d_reco_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_truth_reco = HW2D(dummy_selected_truth_reco, univs, clear_bands);
    


    MH2D* dummy_selected_data_reco =
        new MH2D(Form("selected_data2d_reco_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_data_reco = HW2D(dummy_selected_data_reco, univs, clear_bands);

    delete dummy_selected_mc_reco;
    delete dummy_selected_truth_reco;
    delete dummy_selected_data_reco;
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f,bool isMC) const {
    f.cd();
       if(isMC) m_selected_mc_reco.hist->Write();
       else m_selected_truth_reco.hist->Write();
    // selected mc reco
  }
};
}  // namespace Var2DLoop


#endif  // VARIABLE_H
