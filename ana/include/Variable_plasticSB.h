#ifndef VARIABLE_H
#define VARIABLE_H

#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"
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
  HW m_selected_data_reco_US, m_selected_data_reco_DS, m_selected_data_reco_tgt;

  HW m_hists_US_Fe, m_hists_US_Pb, m_hists_US_C, m_hists_US_other, m_hists_US_regUS, m_hists_US_regDS; 
  HW m_hists_DS_Fe, m_hists_DS_Pb, m_hists_DS_C, m_hists_DS_other, m_hists_DS_regUS, m_hists_DS_regDS; 
  HW m_hists_tgt_Fe, m_hists_tgt_Pb, m_hists_tgt_C, m_hists_tgt_other, m_hists_tgt_regUS, m_hists_tgt_regDS; 

  // HISTFOLIO
  // selected mc reco - signal background histfolio
  PlotUtils::HistFolio<MH1D> m_selected_mc_US, m_selected_mc_DS, m_selected_mc_tgt ;
  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHistograms(T univs, int targetID, int targetZ) {
  //void InitializeAllHistograms(T univs) {
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    const bool clear_bands = true;  // we want empty histograms

    // HISTWRAPPER
    // selected mc reco histwrapper
 
    MH1D* selected_data_reco_US = new MH1D(Form("selected_data_reco_US_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco_US = HW(selected_data_reco_US, univs, clear_bands);

    MH1D* selected_data_reco_DS = new MH1D(Form("selected_data_reco_DS_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco_DS = HW(selected_data_reco_DS, univs, clear_bands);

    MH1D* selected_data_reco_tgt = new MH1D(Form("selected_data_reco_tgt_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco_tgt = HW(selected_data_reco_tgt, univs, clear_bands);

  // HISTFOLIO
  // selected mc reco - signal background histfolio
  std::cout<<"TargetID in varibale.h = "<<targetID<<std::endl;
  std::cout<<"TargetZ in varibale.h = "<<targetZ<<std::endl;

  string trueZ;
  if( targetZ == 26 ) trueZ = "Iron";
  if( targetZ == 82 ) trueZ = "Lead";
  if( targetZ == 6 ) trueZ = "Carbon";
  std::cout<<"I am going to make plots for "<<trueZ<<std::endl;
//-------------------------------------------------------------------------------------------------------------------------
{
    MH1D* hists_US_Fe = new MH1D(Form("US_Fe_%s_%s", trueZ.c_str(), name), "Iron", GetNBins(), bins.data());
    m_hists_US_Fe = HW(hists_US_Fe, univs, clear_bands);

    MH1D* hists_US_Pb = new MH1D(Form("US_Pb_%s_%s", trueZ.c_str(), name), "Lead", GetNBins(), bins.data());
    m_hists_US_Pb = HW(hists_US_Pb, univs, clear_bands);

    MH1D* hists_US_C = new MH1D(Form("US_C_%s_%s", trueZ.c_str(), name), "Carbon", GetNBins(), bins.data());
    m_hists_US_C = HW(hists_US_C, univs, clear_bands);

    MH1D* hists_US_other = new MH1D(Form("US_other_%s_%s", trueZ.c_str(), name), "Other Target", GetNBins(), bins.data());
    m_hists_US_other = HW(hists_US_other, univs, clear_bands);

    MH1D* hists_US_regUS = new MH1D(Form("US_regUS_%s_%s", trueZ.c_str(), name), "Upstream Scintillator", GetNBins(), bins.data());
    m_hists_US_regUS = HW(hists_US_regUS, univs, clear_bands);

    MH1D* hists_US_regDS = new MH1D(Form("US_regDS_%s_%s", trueZ.c_str(), name), "Downstream Scintillator", GetNBins(), bins.data());
    m_hists_US_regDS = HW(hists_US_regDS, univs, clear_bands);

    MH1D* hists_DS_Fe = new MH1D(Form("DS_Fe_%s_%s", trueZ.c_str(), name), "Iron", GetNBins(), bins.data());
    m_hists_DS_Fe = HW(hists_DS_Fe, univs, clear_bands);

    MH1D* hists_DS_Pb = new MH1D(Form("DS_Pb_%s_%s", trueZ.c_str(), name), "Lead", GetNBins(), bins.data());
    m_hists_DS_Pb = HW(hists_DS_Pb, univs, clear_bands);

    MH1D* hists_DS_C = new MH1D(Form("DS_C_%s_%s", trueZ.c_str(), name), "Carbon", GetNBins(), bins.data());
    m_hists_DS_C = HW(hists_DS_C, univs, clear_bands);

    MH1D* hists_DS_other = new MH1D(Form("DS_other_%s_%s", trueZ.c_str(), name), "Other Target", GetNBins(), bins.data());
    m_hists_DS_other = HW(hists_DS_other, univs, clear_bands);

    MH1D* hists_DS_regUS = new MH1D(Form("DS_regUS_%s_%s", trueZ.c_str(), name), "Upstream Scintillator", GetNBins(), bins.data());
    m_hists_DS_regUS = HW(hists_DS_regUS, univs, clear_bands);

    MH1D* hists_DS_regDS = new MH1D(Form("DS_regDS_%s_%s", trueZ.c_str(), name), "Downstream Scintillator", GetNBins(), bins.data());
    m_hists_DS_regDS = HW(hists_DS_regDS, univs, clear_bands);

    MH1D* hists_tgt_Fe = new MH1D(Form("Tgt_Fe_%s_%s", trueZ.c_str(), name), "Iron", GetNBins(), bins.data());
    m_hists_tgt_Fe = HW(hists_tgt_Fe, univs, clear_bands);

    MH1D* hists_tgt_Pb = new MH1D(Form("Tgt_Pb_%s_%s", trueZ.c_str(), name), "Lead", GetNBins(), bins.data());
    m_hists_tgt_Pb = HW(hists_tgt_Pb, univs, clear_bands);

    MH1D* hists_tgt_C = new MH1D(Form("Tgt_C_%s_%s", trueZ.c_str(), name), "Carbon", GetNBins(), bins.data());
    m_hists_tgt_C = HW(hists_tgt_C, univs, clear_bands);

    MH1D* hists_tgt_other = new MH1D(Form("Tgt_other_%s_%s", trueZ.c_str(), name), "Other Target", GetNBins(), bins.data());
    m_hists_tgt_other = HW(hists_tgt_other, univs, clear_bands);

    MH1D* hists_tgt_regUS = new MH1D(Form("Tgt_regUS_%s_%s", trueZ.c_str(), name), "Upstream Scintillator", GetNBins(), bins.data());
    m_hists_tgt_regUS = HW(hists_tgt_regUS, univs, clear_bands);

    MH1D* hists_tgt_regDS = new MH1D(Form("Tgt_regDS_%s_%s", trueZ.c_str(), name), "Downstream Scintillator", GetNBins(), bins.data());
    m_hists_tgt_regDS = HW(hists_tgt_regDS, univs, clear_bands);
}
//-------------------------------------------------------------------------------------------------------------------------
delete selected_data_reco_US;
delete selected_data_reco_DS;
delete selected_data_reco_tgt;

  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
 void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
if(isMC){
  m_hists_US_Fe.hist->Write();
  m_hists_US_Pb.hist->Write();
  m_hists_US_C.hist->Write();
  m_hists_US_other.hist->Write();
  m_hists_US_regUS.hist->Write();
  m_hists_US_regDS.hist->Write();
  m_hists_DS_Fe.hist->Write();
  m_hists_DS_Pb.hist->Write();
  m_hists_DS_C.hist->Write();
  m_hists_DS_other.hist->Write();
  m_hists_DS_regUS.hist->Write();
  m_hists_DS_regDS.hist->Write();
  m_hists_tgt_Fe.hist->Write();
  m_hists_tgt_Pb.hist->Write();
  m_hists_tgt_C.hist->Write();
  m_hists_tgt_other.hist->Write();
  m_hists_tgt_regUS.hist->Write();
  m_hists_tgt_regDS.hist->Write();
    }
//    else m_selected_data_reco.hist->Write();

    // selected mc  histfolio fir Hist Stacking
  // if(isMC){
 //m_selected_mc_US.WriteToFile(f);
 //m_selected_mc_DS.WriteToFile(f);
 //m_selected_mc_tgt.WriteToFile(f);
  //}
   else {
  m_selected_data_reco_US.hist->Write();
  m_selected_data_reco_DS.hist->Write();
  m_selected_data_reco_tgt.hist->Write();
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
//  HW2D m_selected_mc_reco,m_selected_data_reco;
  HW2D m_selected_data_reco_US, m_selected_data_reco_DS, m_selected_data_reco_tgt;

  HW2D m_hists_US_Fe, m_hists_US_Pb, m_hists_US_C, m_hists_US_other, m_hists_US_regUS, m_hists_US_regDS; 
  HW2D m_hists_DS_Fe, m_hists_DS_Pb, m_hists_DS_C, m_hists_DS_other, m_hists_DS_regUS, m_hists_DS_regDS; 
  HW2D m_hists_tgt_Fe, m_hists_tgt_Pb, m_hists_tgt_C, m_hists_tgt_other, m_hists_tgt_regUS, m_hists_tgt_regDS; 


  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHistograms(T univs) {
    const bool clear_bands = true;  // we want empty histograms
    const char* name = GetName().c_str();
    
    MH2D* dummy_selected_data_reco_US =
        new MH2D(Form("selected_data_reco_US_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_data_reco_US = HW2D(dummy_selected_data_reco_US, univs, clear_bands);
  
    MH2D* dummy_selected_data_reco_DS =
        new MH2D(Form("selected_data_reco_DS_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_data_reco_DS = HW2D(dummy_selected_data_reco_DS, univs, clear_bands);

    MH2D* dummy_selected_data_reco_tgt =
        new MH2D(Form("selected_data_reco_tgt_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_data_reco_tgt = HW2D(dummy_selected_data_reco_tgt, univs, clear_bands);

    MH2D* dummy_selected_mc_US_Fe =
        new MH2D(Form("selected_mc_US_Fe_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_US_Fe = HW2D(dummy_selected_mc_US_Fe, univs, clear_bands);
  
    MH2D* dummy_selected_mc_US_Pb =
        new MH2D(Form("selected_mc_US_Pb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_US_Pb = HW2D(dummy_selected_mc_US_Pb, univs, clear_bands);

    MH2D* dummy_selected_mc_US_C =
        new MH2D(Form("selected_mc_US_C_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_US_C = HW2D(dummy_selected_mc_US_C, univs, clear_bands);

    MH2D* dummy_selected_mc_US_other =
        new MH2D(Form("selected_mc_US_other_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_US_other = HW2D(dummy_selected_mc_US_other, univs, clear_bands);

    MH2D* dummy_selected_mc_US_regUS =
        new MH2D(Form("selected_mc_US_regUS_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_US_regUS = HW2D(dummy_selected_mc_US_regUS, univs, clear_bands);

    MH2D* dummy_selected_mc_US_regDS =
        new MH2D(Form("selected_mc_US_regDS_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_US_regDS = HW2D(dummy_selected_mc_US_regDS, univs, clear_bands);

    MH2D* dummy_selected_mc_DS_Fe =
        new MH2D(Form("selected_mc_DS_Fe_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_DS_Fe = HW2D(dummy_selected_mc_DS_Fe, univs, clear_bands);
  
    MH2D* dummy_selected_mc_DS_Pb =
        new MH2D(Form("selected_mc_DS_Pb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_DS_Pb = HW2D(dummy_selected_mc_DS_Pb, univs, clear_bands);

    MH2D* dummy_selected_mc_DS_C =
        new MH2D(Form("selected_mc_DS_C_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_DS_C = HW2D(dummy_selected_mc_DS_C, univs, clear_bands);

    MH2D* dummy_selected_mc_DS_other =
        new MH2D(Form("selected_mc_DS_other_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_DS_other = HW2D(dummy_selected_mc_DS_other, univs, clear_bands);

    MH2D* dummy_selected_mc_DS_regUS =
        new MH2D(Form("selected_mc_DS_regUS_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_DS_regUS = HW2D(dummy_selected_mc_DS_regUS, univs, clear_bands);

    MH2D* dummy_selected_mc_DS_regDS =
        new MH2D(Form("selected_mc_DS_regDS_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_DS_regDS = HW2D(dummy_selected_mc_DS_regDS, univs, clear_bands);

    MH2D* dummy_selected_mc_tgt_Fe =
        new MH2D(Form("selected_mc_tgt_Fe_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_tgt_Fe = HW2D(dummy_selected_mc_tgt_Fe, univs, clear_bands);
  
    MH2D* dummy_selected_mc_tgt_Pb =
        new MH2D(Form("selected_mc_tgt_Pb_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_tgt_Pb = HW2D(dummy_selected_mc_tgt_Pb, univs, clear_bands);

    MH2D* dummy_selected_mc_tgt_C =
        new MH2D(Form("selected_mc_tgt_C_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_tgt_C = HW2D(dummy_selected_mc_tgt_C, univs, clear_bands);

    MH2D* dummy_selected_mc_tgt_other =
        new MH2D(Form("selected_mc_tgt_other_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_tgt_other = HW2D(dummy_selected_mc_tgt_other, univs, clear_bands);

    MH2D* dummy_selected_mc_tgt_regUS =
        new MH2D(Form("selected_mc_tgt_regUS_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_tgt_regUS = HW2D(dummy_selected_mc_tgt_regUS, univs, clear_bands);

    MH2D* dummy_selected_mc_tgt_regDS =
        new MH2D(Form("selected_mc_tgt_regDS_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_hists_tgt_regDS = HW2D(dummy_selected_mc_tgt_regDS, univs, clear_bands);


    delete dummy_selected_mc_US_Fe;
    delete dummy_selected_mc_US_Pb;
    delete dummy_selected_mc_US_C;
    delete dummy_selected_mc_US_other;
    delete dummy_selected_mc_US_regUS;
    delete dummy_selected_mc_US_regDS;
    delete dummy_selected_mc_DS_Fe;
    delete dummy_selected_mc_DS_Pb;
    delete dummy_selected_mc_DS_C;
    delete dummy_selected_mc_DS_other;
    delete dummy_selected_mc_DS_regUS;
    delete dummy_selected_mc_DS_regDS;
    delete dummy_selected_mc_tgt_Fe;
    delete dummy_selected_mc_tgt_Pb;
    delete dummy_selected_mc_tgt_C;
    delete dummy_selected_mc_tgt_other;
    delete dummy_selected_mc_tgt_regUS;
    delete dummy_selected_mc_tgt_regDS;

  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ 
           //m_selected_mc_reco.hist->Write();
           //m_selected_mc_lowW.hist->Write();
  m_hists_US_Fe.hist->Write();
  m_hists_US_Pb.hist->Write();
  m_hists_US_C.hist->Write();
  m_hists_US_other.hist->Write();
  m_hists_US_regUS.hist->Write();
  m_hists_US_regDS.hist->Write();
  m_hists_DS_Fe.hist->Write();
  m_hists_DS_Pb.hist->Write();
  m_hists_DS_C.hist->Write();
  m_hists_DS_other.hist->Write();
  m_hists_DS_regUS.hist->Write();
  m_hists_DS_regDS.hist->Write();
  m_hists_tgt_Fe.hist->Write();
  m_hists_tgt_Pb.hist->Write();
  m_hists_tgt_C.hist->Write();
  m_hists_tgt_other.hist->Write();
  m_hists_tgt_regUS.hist->Write();
  m_hists_tgt_regDS.hist->Write();
        }
        else {   
  m_selected_data_reco_US.hist->Write();
  m_selected_data_reco_DS.hist->Write();
  m_selected_data_reco_tgt.hist->Write();
        }
//       else m_selected_data_reco.hist->Write();
    // selected mc reco
  }
};
}  // namespace Var2DLoop




#endif  // VARIABLE_H
