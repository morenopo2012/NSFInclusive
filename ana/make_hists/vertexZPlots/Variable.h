#ifndef VARIABLE_H
#define VARIABLE_H

#include "../../../NUKECCSRC/ana_common/include/CVUniverse.h"
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
    HW m_selected_mc_reco_water;
    HW m_selected_mc_reco_carbon;
    HW m_selected_mc_reco_iron;
    HW m_selected_mc_reco_lead;
    HW m_selected_mc_reco_scintillator;
    HW m_selected_data_reco;

    HW m_selected_data_reco_water;
    HW m_selected_data_reco_carbon;
    HW m_selected_data_reco_iron;
    HW m_selected_data_reco_lead;
    HW m_selected_data_reco_scintillator;
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
      MH1D* dummy_selected_mc_reco_water = new MH1D(Form("selected_mc_reco_water_%s", name), name, GetNBins(), bins.data());
      m_selected_mc_reco_water = HW(dummy_selected_mc_reco_water, univs, clear_bands);

      MH1D* dummy_selected_mc_reco_carbon = new MH1D(Form("selected_mc_reco_carbon_%s", name), name, GetNBins(), bins.data());
      m_selected_mc_reco_carbon = HW(dummy_selected_mc_reco_carbon, univs, clear_bands);

      MH1D* dummy_selected_mc_reco_iron = new MH1D(Form("selected_mc_reco_iron_%s", name), name, GetNBins(), bins.data());
      m_selected_mc_reco_iron = HW(dummy_selected_mc_reco_iron, univs, clear_bands);

      MH1D* dummy_selected_mc_reco_lead = new MH1D(Form("selected_mc_reco_lead_%s", name), name, GetNBins(), bins.data());
      m_selected_mc_reco_lead = HW(dummy_selected_mc_reco_lead, univs, clear_bands);

      MH1D* dummy_selected_mc_reco_scintillator = new MH1D(Form("selected_mc_reco_scintillator_%s", name), name, GetNBins(), bins.data());
      m_selected_mc_reco_scintillator = HW(dummy_selected_mc_reco_scintillator, univs, clear_bands);



      MH1D* dummy_selected_data_reco_water = new MH1D(Form("selected_data_reco_water_%s", name), name, GetNBins(), bins.data());
      m_selected_data_reco_water = HW(dummy_selected_data_reco_water, univs, clear_bands);

      MH1D* dummy_selected_data_reco_carbon = new MH1D(Form("selected_data_reco_carbon_%s", name), name, GetNBins(), bins.data());
      m_selected_data_reco_carbon = HW(dummy_selected_data_reco_carbon, univs, clear_bands);

      MH1D* dummy_selected_data_reco_iron = new MH1D(Form("selected_data_reco_iron_%s", name), name, GetNBins(), bins.data());
      m_selected_data_reco_iron = HW(dummy_selected_data_reco_iron, univs, clear_bands);

      MH1D* dummy_selected_data_reco_lead = new MH1D(Form("selected_data_reco_lead_%s", name), name, GetNBins(), bins.data());
      m_selected_data_reco_lead = HW(dummy_selected_data_reco_lead, univs, clear_bands);

      MH1D* dummy_selected_data_reco_scintillator = new MH1D(Form("selected_data_reco_scintillator_%s", name), name, GetNBins(), bins.data());
      m_selected_data_reco_scintillator = HW(dummy_selected_data_reco_scintillator, univs, clear_bands);

      MH1D* dummy_selected_data_reco = new MH1D(Form("selected_data_reco_%s", name), name, GetNBins(), bins.data());
      m_selected_data_reco = HW(dummy_selected_data_reco, univs, clear_bands);
    
      delete dummy_selected_mc_reco_water;
      delete dummy_selected_mc_reco_carbon;
      delete dummy_selected_mc_reco_iron;
      delete dummy_selected_mc_reco_lead;
      delete dummy_selected_mc_reco_scintillator;

      delete dummy_selected_data_reco_water;
      delete dummy_selected_data_reco_carbon;
      delete dummy_selected_data_reco_iron;
      delete dummy_selected_data_reco_lead;
      delete dummy_selected_data_reco_scintillator;
    }

    //=======================================================================================
    // WRITE ALL HISTOGRAMS
    //=======================================================================================
    void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
      f.cd();

      // selected mc reco
      if(isMC) {
        m_selected_mc_reco_water.hist->Write();
        m_selected_mc_reco_carbon.hist->Write();
        m_selected_mc_reco_iron.hist->Write();
        m_selected_mc_reco_lead.hist->Write();
        m_selected_mc_reco_scintillator.hist->Write();
      }
    else{
      m_selected_data_reco.hist->Write();
      m_selected_data_reco_water.hist->Write();
      m_selected_data_reco_carbon.hist->Write();
      m_selected_data_reco_iron.hist->Write();
      m_selected_data_reco_lead.hist->Write();
      m_selected_data_reco_scintillator.hist->Write();
    }
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
  //HW2D m_selected_mc_reco,m_selected_data_reco;

  HW2D m_selected_mc_truth_trackerC, m_selected_mc_truth_t25fe, m_selected_mc_truth_t3c, m_selected_mc_truth_t15fe;
  /*
  HW2Dm_selected_mc_truth_waterO, m_selected_mc_truth_t1fe, m_selected_mc_truth_t1pb;
  HW2D m_selected_mc_truth_t2fe, m_selected_mc_truth_t2pb, m_selected_mc_truth_t3fe, m_selected_mc_truth_t3pb, m_selected_mc_truth_t3c;
  HW2D m_selected_mc_truth_t4pb, m_selected_mc_truth_t5fe, m_selected_mc_truth_t5pb;
  HW2D m_selected_mc_truth_t25fe, m_selected_mc_truth_t25pb;
  HW2D m_selected_mc_truth_t15fe, m_selected_mc_truth_t15pb;
  */
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
    

    MH2D* dummy_selected_mc_truth_trackerC = new MH2D(Form("selected_mc_truth_trackerC_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_trackerC = HW2D(dummy_selected_mc_truth_trackerC, univs, clear_bands);
    /*
    MH2D* dummy_selected_mc_truth_waterO = new MH2D(Form("selected_mc_truth_waterO_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_waterO = HW2D(dummy_selected_mc_truth_waterO, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_t1fe = new MH2D(Form("selected_mc_truth_t1fe_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t1fe = HW2D(dummy_selected_mc_truth_t1fe, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_t1pb = new MH2D(Form("selected_mc_truth_t1pb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t1pb = HW2D(dummy_selected_mc_truth_t1pb, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_t2fe = new MH2D(Form("selected_mc_truth_t2fe_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t2fe = HW2D(dummy_selected_mc_truth_t2fe, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_t2pb = new MH2D(Form("selected_mc_truth_t2pb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t2pb = HW2D(dummy_selected_mc_truth_t2pb, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_t3fe = new MH2D(Form("selected_mc_truth_t3fe_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t3fe = HW2D(dummy_selected_mc_truth_t3fe, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_t3pb = new MH2D(Form("selected_mc_truth_t3pb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t3pb = HW2D(dummy_selected_mc_truth_t3pb, univs, clear_bands);
    */
    MH2D* dummy_selected_mc_truth_t3c = new MH2D(Form("selected_mc_truth_t3c_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t3c = HW2D(dummy_selected_mc_truth_t3c, univs, clear_bands);
    /*
    MH2D* dummy_selected_mc_truth_t4pb = new MH2D(Form("selected_mc_truth_t4pb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t4pb = HW2D(dummy_selected_mc_truth_t4pb, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_t5fe = new MH2D(Form("selected_mc_truth_t5fe_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t5fe = HW2D(dummy_selected_mc_truth_t5fe, univs, clear_bands);

    MH2D* dummy_selected_mc_truth_t5pb = new MH2D(Form("selected_mc_truth_t5pb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t5pb = HW2D(dummy_selected_mc_truth_t5pb, univs, clear_bands);
*/
    MH2D* dummy_selected_mc_truth_t25fe = new MH2D(Form("selected_mc_truth_t25fe_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t25fe = HW2D(dummy_selected_mc_truth_t25fe, univs, clear_bands);
 /*
    MH2D* dummy_selected_mc_truth_t25pb = new MH2D(Form("selected_mc_truth_t25pb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t25pb = HW2D(dummy_selected_mc_truth_t25pb, univs, clear_bands);
*/
    MH2D* dummy_selected_mc_truth_t15fe = new MH2D(Form("selected_mc_truth_t15fe_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t15fe = HW2D(dummy_selected_mc_truth_t15fe, univs, clear_bands);
 /*
    MH2D* dummy_selected_mc_truth_t15pb = new MH2D(Form("selected_mc_truth_t15pb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truth_t15pb = HW2D(dummy_selected_mc_truth_t15pb, univs, clear_bands);
    */
    delete dummy_selected_mc_truth_trackerC;
    delete dummy_selected_mc_truth_t3c;
    delete dummy_selected_mc_truth_t25fe;
    delete dummy_selected_mc_truth_t15fe;
    /*
    delete dummy_selected_mc_truth_waterO;
    delete dummy_selected_mc_truth_t1fe;
    delete dummy_selected_mc_truth_t1pb;
    delete dummy_selected_mc_truth_t2fe;
    delete dummy_selected_mc_truth_t2pb;
    delete dummy_selected_mc_truth_t3fe;
    delete dummy_selected_mc_truth_t3pb;
    delete dummy_selected_mc_truth_t3c;
    delete dummy_selected_mc_truth_t4pb;
    delete dummy_selected_mc_truth_t5fe;
    delete dummy_selected_mc_truth_t5pb;
    delete dummy_selected_mc_truth_t25pb;
    delete dummy_selected_mc_truth_t15fe;
    delete dummy_selected_mc_truth_t15pb;
    */
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f,bool isMC) const {
    f.cd();
       if(isMC) {
        m_selected_mc_truth_trackerC.hist->Write();
        m_selected_mc_truth_t3c.hist->Write();
        m_selected_mc_truth_t25fe.hist->Write();
        m_selected_mc_truth_t15fe.hist->Write();
        /*
        m_selected_mc_truth_waterO.hist->Write();
        m_selected_mc_truth_t1fe.hist->Write();
        m_selected_mc_truth_t1pb.hist->Write();
        m_selected_mc_truth_t2fe.hist->Write();
        m_selected_mc_truth_t2pb.hist->Write();
        m_selected_mc_truth_t3fe.hist->Write();
        m_selected_mc_truth_t3pb.hist->Write();
        m_selected_mc_truth_t3c.hist->Write();
        m_selected_mc_truth_t4pb.hist->Write();
        m_selected_mc_truth_t5fe.hist->Write();
        m_selected_mc_truth_t5pb.hist->Write();
        m_selected_mc_truth_t25fe.hist->Write();
        m_selected_mc_truth_t25pb.hist->Write();
        m_selected_mc_truth_t15pb.hist->Write();
        */
      }
    // selected mc reco
  }
};
}  // namespace Var2DLoop


#endif  // VARIABLE_H
