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
  HW m_selected_data_reco_tgt, m_selected_data_reco_trans, m_selected_data_reco_contin;

  // HISTFOLIO
  // selected mc reco - signal background histfolio
  PlotUtils::HistFolio<MH1D> m_selected_mc_trans, m_selected_mc_contin, m_selected_mc_sig ;
  HW m_hists_trans_mc_CH_DS, m_hists_contin_mc_CH_DS;
  HW m_hists_trans_mc_CH_US, m_hists_contin_mc_CH_US;
  HW m_hists_trans_in_trans, m_hists_contin_in_trans, m_hists_signal_in_trans;
  HW m_hists_trans_in_contin, m_hists_contin_in_contin, m_hists_signal_in_contin;
  HW m_hists_true_signal;

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
    //MH1D* dummy_selected_mc_reco = new MH1D(Form("selected_mc_reco_%s", name), name, GetNBins(), bins.data());
   // m_selected_mc_reco = HW(dummy_selected_mc_reco, univs, clear_bands);

//    MH1D* dummy_selected_mc_reco_US = new MH1D(Form("selected_mc_reco_US_%s", name), name, GetNBins(), bins.data());
//    m_selected_mc_reco_US = HW(dummy_selected_mc_reco_US, univs, clear_bands);

//    MH1D* dummy_selected_mc_reco_DS = new MH1D(Form("selected_mc_reco_DS_%s", name), name, GetNBins(), bins.data());
//    m_selected_mc_reco_DS = HW(dummy_selected_mc_reco_DS, univs, clear_bands);

//    MH1D* dummy_selected_data_reco = new MH1D(Form("selected_data_reco_%s", name), name, GetNBins(), bins.data());
//    m_selected_data_reco = HW(dummy_selected_data_reco, univs, clear_bands);

//    MH1D* selected_data_reco_sb = new MH1D(Form("selected_data_reco_sb_%s", name), name, GetNBins(), bins.data());
//    m_selected_data_reco_sb = HW(selected_data_reco_sb, univs, clear_bands);
 
    //MH1D* selected_data_reco_US = new MH1D(Form("selected_data_reco_US_%s", name), name, GetNBins(), bins.data());
    //m_selected_data_reco_US = HW(selected_data_reco_US, univs, clear_bands);

    //MH1D* selected_data_reco_DS = new MH1D(Form("selected_data_reco_DS_%s", name), name, GetNBins(), bins.data());
    //m_selected_data_reco_DS = HW(selected_data_reco_DS, univs, clear_bands);

    MH1D* selected_data_reco_tgt = new MH1D(Form("selected_data_reco_tgt_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco_tgt = HW(selected_data_reco_tgt, univs, clear_bands);

    MH1D* selected_data_reco_trans = new MH1D(Form("selected_data_reco_trans_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco_trans = HW(selected_data_reco_trans, univs, clear_bands);

    MH1D* selected_data_reco_contin = new MH1D(Form("selected_data_reco_contin_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco_contin = HW(selected_data_reco_contin, univs, clear_bands);

    MH1D* hists_trans_mc_CH_US = new MH1D(Form("hists_trans_mc_CH_US_%s", name), name, GetNBins(), bins.data());
    m_hists_trans_mc_CH_US = HW(hists_trans_mc_CH_US, univs, clear_bands);

    MH1D* hists_contin_mc_CH_US = new MH1D(Form("hists_contin_mc_CH_US_%s", name), name, GetNBins(), bins.data());
    m_hists_contin_mc_CH_US = HW(hists_contin_mc_CH_US, univs, clear_bands);

    MH1D* hists_trans_mc_CH_DS = new MH1D(Form("hists_trans_mc_CH_DS_%s", name), name, GetNBins(), bins.data());
    m_hists_trans_mc_CH_DS = HW(hists_trans_mc_CH_DS, univs, clear_bands);

    MH1D* hists_contin_mc_CH_DS = new MH1D(Form("hists_contin_mc_CH_DS_%s", name), name, GetNBins(), bins.data());
    m_hists_contin_mc_CH_DS = HW(hists_contin_mc_CH_DS, univs, clear_bands);

    MH1D* hists_trans_in_trans = new MH1D(Form("hists_trans_in_trans_%s", name), name, GetNBins(), bins.data());
    m_hists_trans_in_trans = HW(hists_trans_in_trans, univs, clear_bands);

    MH1D* hists_contin_in_trans = new MH1D(Form("hists_contin_in_trans_%s", name), name, GetNBins(), bins.data());
    m_hists_contin_in_trans = HW(hists_contin_in_trans, univs, clear_bands);

    MH1D* hists_signal_in_trans = new MH1D(Form("hists_signal_in_trans_%s", name), name, GetNBins(), bins.data());
    m_hists_signal_in_trans = HW(hists_signal_in_trans, univs, clear_bands);

    MH1D* hists_trans_in_contin = new MH1D(Form("hists_trans_in_contin_%s", name), name, GetNBins(), bins.data());
    m_hists_trans_in_contin = HW(hists_trans_in_contin, univs, clear_bands);

    MH1D* hists_contin_in_contin = new MH1D(Form("hists_contin_in_contin_%s", name), name, GetNBins(), bins.data());
    m_hists_contin_in_contin = HW(hists_contin_in_contin, univs, clear_bands);

    MH1D* hists_signal_in_contin = new MH1D(Form("hists_signal_in_contin_%s", name), name, GetNBins(), bins.data());
    m_hists_signal_in_contin = HW(hists_signal_in_contin, univs, clear_bands);

    MH1D* hists_true_signal = new MH1D(Form("hists_true_signal_%s", name), name, GetNBins(), bins.data());
    m_hists_true_signal = HW(hists_true_signal, univs, clear_bands);


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
/*
{
    MH1D* hists_US_Fe = new MH1D(Form("US_Fe_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_US_Fe = HW(hists_US_Fe, univs, clear_bands);

    MH1D* hists_US_Pb = new MH1D(Form("US_Pb_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_US_Pb = HW(hists_US_Pb, univs, clear_bands);

    MH1D* hists_US_C = new MH1D(Form("US_C_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_US_C = HW(hists_US_C, univs, clear_bands);

    MH1D* hists_US_other = new MH1D(Form("US_other_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_US_other = HW(hists_US_other, univs, clear_bands);

    MH1D* hists_US_regUS = new MH1D(Form("US_regUS_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_US_regUS = HW(hists_US_regUS, univs, clear_bands);

    MH1D* hists_US_regDS = new MH1D(Form("US_regDS_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_US_regDS = HW(hists_US_regDS, univs, clear_bands);

    MH1D* hists_DS_Fe = new MH1D(Form("DS_Fe_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_DS_Fe = HW(hists_DS_Fe, univs, clear_bands);

    MH1D* hists_DS_Pb = new MH1D(Form("DS_Pb_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_DS_Pb = HW(hists_DS_Pb, univs, clear_bands);

    MH1D* hists_DS_C = new MH1D(Form("DS_C_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_DS_C = HW(hists_DS_C, univs, clear_bands);

    MH1D* hists_DS_other = new MH1D(Form("DS_other_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_DS_other = HW(hists_DS_other, univs, clear_bands);

    MH1D* hists_DS_regUS = new MH1D(Form("DS_regUS_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_DS_regUS = HW(hists_DS_regUS, univs, clear_bands);

    MH1D* hists_DS_regDS = new MH1D(Form("DS_regDS_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_DS_regDS = HW(hists_DS_regDS, univs, clear_bands);

    MH1D* hists_tgt_Fe = new MH1D(Form("Tgt_Fe_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_tgt_Fe = HW(hists_tgt_Fe, univs, clear_bands);

    MH1D* hists_tgt_Pb = new MH1D(Form("Tgt_Pb_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_tgt_Pb = HW(hists_tgt_Pb, univs, clear_bands);

    MH1D* hists_tgt_C = new MH1D(Form("Tgt_C_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_tgt_C = HW(hists_tgt_C, univs, clear_bands);

    MH1D* hists_tgt_other = new MH1D(Form("Tgt_other_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_tgt_other = HW(hists_tgt_other, univs, clear_bands);

    MH1D* hists_tgt_regUS = new MH1D(Form("Tgt_regUS_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_tgt_regUS = HW(hists_tgt_regUS, univs, clear_bands);

    MH1D* hists_tgt_regDS = new MH1D(Form("Tgt_regDS_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_hists_tgt_regDS = HW(hists_tgt_regDS, univs, clear_bands);
}
*/
//-------------------------------------------------------------------------------------------------------------------------
/*
   m_selected_mc_US = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("US_t%d_z%d_%s", targetID, targetZ, name), name, GetNBins(), bins.data());
    if(targetZ == 26)m_selected_mc_US.AddComponentHist("Fe");
    if(targetZ == 82)m_selected_mc_US.AddComponentHist("Pb");
    if(targetZ == 6)m_selected_mc_US.AddComponentHist("C");
    m_selected_mc_US.AddComponentHist("Other");
    m_selected_mc_US.AddComponentHist("regUS");
    m_selected_mc_US.AddComponentHist("regDS");
    
   m_selected_mc_DS = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("DS_t%d_z%d_%s", targetID, targetZ, name), name, GetNBins(), bins.data());
    if(targetZ == 26)m_selected_mc_DS.AddComponentHist("Fe");
    if(targetZ == 82)m_selected_mc_DS.AddComponentHist("Pb");
    if(targetZ == 6)m_selected_mc_DS.AddComponentHist("C");
    m_selected_mc_DS.AddComponentHist("Other");
    m_selected_mc_DS.AddComponentHist("regUS");
    m_selected_mc_DS.AddComponentHist("regDS");
 
   m_selected_mc_tgt = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("Target_t%d_z%d_%s", targetID, targetZ, name), name, GetNBins(), bins.data());
    if(targetZ == 26)m_selected_mc_tgt.AddComponentHist("Fe");
    if(targetZ == 82)m_selected_mc_tgt.AddComponentHist("Pb");
    if(targetZ == 6)m_selected_mc_tgt.AddComponentHist("C");
    m_selected_mc_tgt.AddComponentHist("Other");
    m_selected_mc_tgt.AddComponentHist("regUS");
    m_selected_mc_tgt.AddComponentHist("regDS");
*/ 
/*   m_selected_mc_trans = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("Trans_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_selected_mc_trans.AddComponentHist("True_Trans");
    m_selected_mc_trans.AddComponentHist("True_Contin");
    m_selected_mc_trans.AddComponentHist("True_DIS");

   m_selected_mc_contin = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("Contin_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_selected_mc_contin.AddComponentHist("True_Trans");
    m_selected_mc_contin.AddComponentHist("True_Contin");
    m_selected_mc_contin.AddComponentHist("True_DIS");


   m_selected_mc_sig = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("Signal_%s_%s", trueZ.c_str(), name), name, GetNBins(), bins.data());
    m_selected_mc_sig.AddComponentHist("Signal");

*/
delete selected_data_reco_tgt; 
delete selected_data_reco_trans;
delete selected_data_reco_contin;
delete hists_trans_mc_CH_US;
delete hists_contin_mc_CH_US; 
delete hists_trans_mc_CH_DS;
delete hists_contin_mc_CH_DS;
delete hists_trans_in_trans;
delete hists_contin_in_trans;
delete hists_signal_in_trans;
delete hists_trans_in_contin;
delete hists_contin_in_contin;
delete hists_signal_in_contin;
delete hists_true_signal;
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
 void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
if(isMC){
//  m_selected_mc_reco.hist->Write();
//  m_selected_mc_reco_US.hist->Write();
//  m_selected_mc_reco_DS.hist->Write();
/*  m_hists_US_Fe.hist->Write();
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
*/
  m_hists_trans_mc_CH_US.hist->Write();
  m_hists_contin_mc_CH_US.hist->Write();
  m_hists_trans_mc_CH_DS.hist->Write();
  m_hists_contin_mc_CH_DS.hist->Write();
  m_hists_true_signal.hist->Write();
  m_hists_trans_in_trans.hist->Write();
  m_hists_contin_in_trans.hist->Write();
  m_hists_signal_in_trans.hist->Write();
  m_hists_trans_in_contin.hist->Write();
  m_hists_contin_in_contin.hist->Write();
  m_hists_signal_in_contin.hist->Write();
    }
//    else m_selected_data_reco.hist->Write();

    // selected mc  histfolio fir Hist Stacking
   //if(isMC){
 //m_selected_mc_sb.WriteToFile(f);
 //m_selected_mc_US.WriteToFile(f);
 //m_selected_mc_DS.WriteToFile(f);
 //m_selected_mc_tgt.WriteToFile(f);
 //m_selected_mc_trans.WriteToFile(f);
 //m_selected_mc_contin.WriteToFile(f);
 //m_selected_mc_sig.WriteToFile(f);
  //}
   else {
  //m_selected_data_reco_sb.hist->Write();
  //m_selected_data_reco_US.hist->Write();
  //m_selected_data_reco_DS.hist->Write();
  m_selected_data_reco_tgt.hist->Write();
  m_selected_data_reco_trans.hist->Write();
  m_selected_data_reco_contin.hist->Write();
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
  HW2D m_selected_mc_lowW, m_selected_mc_lowQ2, m_selected_mc_dis, m_selected_mc_contin_in_trans, m_selected_mc_trans_in_trans, m_selected_mc_signal_in_trans, m_selected_true_signal, m_selected_trans_CH_US, m_selected_contin_CH_US, m_selected_data_reco_trans, m_selected_data_reco_contin;
  HW2D m_selected_mc_q2, m_selected_mc_trueLowW, m_selected_mc_trueQ2, m_selected_mc_truetrans, m_selected_mc_contin_in_contin, m_selected_mc_trans_in_contin, m_selected_mc_signal_in_contin, m_selected_trans_CH_DS, m_selected_contin_CH_DS;

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
//    MH2D* dummy_selected_mc_reco =
//        new MH2D(Form("selected_mc_reco2d_%s", name), name, GetNBinsX(),
//                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
//    m_selected_mc_reco = HW2D(dummy_selected_mc_reco, univs, clear_bands);
/*    
     MH2D* dummy_selected_mc_lowW =
        new MH2D(Form("selected_mc_lowW_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_lowW = HW2D(dummy_selected_mc_lowW, univs, clear_bands);
    
    MH2D* dummy_selected_mc_lowQ2 =
        new MH2D(Form("selected_mc_lowQ2_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_lowQ2 = HW2D(dummy_selected_mc_lowQ2, univs, clear_bands);
 */
    MH2D* dummy_selected_mc_dis =
        new MH2D(Form("selected_mc_dis_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_dis = HW2D(dummy_selected_mc_dis, univs, clear_bands);
  
    MH2D* dummy_selected_mc_q2 =
        new MH2D(Form("selected_mc_q2_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_q2 = HW2D(dummy_selected_mc_q2, univs, clear_bands);
  
    MH2D* dummy_selected_mc_trueLowW =
        new MH2D(Form("selected_mc_trueLowW_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_trueLowW = HW2D(dummy_selected_mc_trueLowW, univs, clear_bands);
  
    MH2D* dummy_selected_mc_trueQ2 =
        new MH2D(Form("selected_mc_trueQ2_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_trueQ2 = HW2D(dummy_selected_mc_trueQ2, univs, clear_bands);
  
    MH2D* dummy_selected_mc_truetrans =
        new MH2D(Form("selected_mc_truetrans_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_truetrans = HW2D(dummy_selected_mc_truetrans, univs, clear_bands);

////////////////*********2D physics sideband stuff starts********/////////
    MH2D* dummy_selected_mc_contin_in_trans =
        new MH2D(Form("selected_mc_contin_in_trans_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_contin_in_trans = HW2D(dummy_selected_mc_contin_in_trans, univs, clear_bands);

    MH2D* dummy_selected_mc_trans_in_trans =
        new MH2D(Form("selected_mc_trans_in_trans_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_trans_in_trans = HW2D(dummy_selected_mc_trans_in_trans, univs, clear_bands);

    MH2D* dummy_selected_mc_signal_in_trans =
        new MH2D(Form("selected_mc_signal_in_trans_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_signal_in_trans = HW2D(dummy_selected_mc_signal_in_trans, univs, clear_bands);

    MH2D* dummy_selected_mc_contin_in_contin =
        new MH2D(Form("selected_mc_contin_in_contin_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_contin_in_contin = HW2D(dummy_selected_mc_contin_in_contin, univs, clear_bands);

    MH2D* dummy_selected_mc_trans_in_contin =
        new MH2D(Form("selected_mc_trans_in_contin_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_trans_in_contin = HW2D(dummy_selected_mc_trans_in_contin, univs, clear_bands);

    MH2D* dummy_selected_mc_signal_in_contin =
        new MH2D(Form("selected_mc_signal_in_contin_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_signal_in_contin = HW2D(dummy_selected_mc_signal_in_contin, univs, clear_bands);


    MH2D* dummy_selected_mc_true_signal =
        new MH2D(Form("selected_true_in_signal_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_true_signal = HW2D(dummy_selected_mc_true_signal, univs, clear_bands);

    MH2D* dummy_selected_mc_trans_CH_US =
        new MH2D(Form("selected_trans_CH_US_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_trans_CH_US = HW2D(dummy_selected_mc_trans_CH_US, univs, clear_bands);

    MH2D* dummy_selected_mc_trans_CH_DS =
        new MH2D(Form("selected_trans_CH_DS_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_trans_CH_DS = HW2D(dummy_selected_mc_trans_CH_DS, univs, clear_bands);

    MH2D* dummy_selected_mc_contin_CH_US =
        new MH2D(Form("selected_contin_CH_US_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_contin_CH_US = HW2D(dummy_selected_mc_contin_CH_US, univs, clear_bands);

    MH2D* dummy_selected_mc_contin_CH_DS =
        new MH2D(Form("selected_contin_CH_DS_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_contin_CH_DS = HW2D(dummy_selected_mc_contin_CH_DS, univs, clear_bands);
    
    MH2D* dummy_selected_data_trans =
        new MH2D(Form("selected_data_trans_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_data_reco_trans = HW2D(dummy_selected_data_trans, univs, clear_bands);

    MH2D* dummy_selected_data_contin =
        new MH2D(Form("selected_data_contin_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_data_reco_contin = HW2D(dummy_selected_data_contin, univs, clear_bands);
////////////////*********2D physics sideband stuff ends********/////////

//    MH2D* dummy_selected_data_reco =
//        new MH2D(Form("selected_data2d_reco_%s", name), name, GetNBinsX(),
//                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
//    m_selected_data_reco = HW2D(dummy_selected_data_reco, univs, clear_bands);

    //delete dummy_selected_mc_reco;
    //delete dummy_selected_mc_lowW;
    //delete dummy_selected_mc_lowQ2;
    delete dummy_selected_mc_dis;
    delete dummy_selected_mc_q2;
    delete dummy_selected_mc_trueLowW;
    delete dummy_selected_mc_trueQ2;
    delete dummy_selected_mc_truetrans;
    delete dummy_selected_mc_trans_in_trans;
    delete dummy_selected_mc_contin_in_trans;
    delete dummy_selected_mc_signal_in_trans;
    delete dummy_selected_mc_trans_in_contin;
    delete dummy_selected_mc_contin_in_contin;
    delete dummy_selected_mc_signal_in_contin;
    delete dummy_selected_mc_true_signal;
    delete dummy_selected_mc_trans_CH_US;
    delete dummy_selected_mc_trans_CH_DS;
    delete dummy_selected_mc_contin_CH_US;
    delete dummy_selected_mc_contin_CH_DS;
    delete dummy_selected_data_trans;
    delete dummy_selected_data_contin; 

   // delete dummy_selected_data_reco;
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ 
           //m_selected_mc_reco.hist->Write();
           //m_selected_mc_lowW.hist->Write();
           //m_selected_mc_lowQ2.hist->Write();
           m_selected_mc_dis.hist->Write();
           m_selected_mc_q2.hist->Write();
           m_selected_mc_trueLowW.hist->Write();
           m_selected_mc_trueQ2.hist->Write();
           m_selected_mc_truetrans.hist->Write();
           m_selected_mc_trans_in_trans.hist->Write();
           m_selected_mc_contin_in_trans.hist->Write();
           m_selected_mc_signal_in_trans.hist->Write();
           m_selected_mc_trans_in_contin.hist->Write();
           m_selected_mc_contin_in_contin.hist->Write();
           m_selected_mc_signal_in_contin.hist->Write();
           m_selected_true_signal.hist->Write();
           m_selected_trans_CH_US.hist->Write();
           m_selected_trans_CH_DS.hist->Write();
           m_selected_contin_CH_US.hist->Write();
           m_selected_contin_CH_DS.hist->Write(); 
        }
        else {   
           m_selected_data_reco_trans.hist->Write();
           m_selected_data_reco_contin.hist->Write(); 
        }
//       else m_selected_data_reco.hist->Write();
    // selected mc reco
  }
};
}  // namespace Var2DLoop


#endif  // VARIABLE_H
