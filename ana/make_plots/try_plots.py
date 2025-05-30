import ROOT
import os,sys
from PlotUtils import MnvH1D,MnvH2D,MnvPlotter
_fnsf = ROOT.TFile("/pnfs/minerva/persistent/users/hhaider/NukeHists/v1/Hists_Energy_MC_t1_z26_Nu_v1_.root","READ")
mc_name = "Enu_CC"
mc_name1 = "Enu_dis"
mc_new = _fnsf.Get(mc_name)
mc_new1 = _fnsf.Get(mc_name1)
canvas = ROOT.TCanvas()
canvas.cd()
mc_new_cv = mc_new.GetCVHistoWithStatError()
mc_new_cv1 = mc_new1.GetCVHistoWithStatError()
#mc_new_cv.SetTitle("NSF")
mc_new_cv.SetLineColor(2)
mc_new_cv.Draw("hist")
mc_new_cv1.Draw("hist same")
#canvas.BuildLegend()
#raw_input()
canvas.Print("trynuE.png") 


mnv = MnvPlotter()
canvas1 = ROOT.TCanvas()
mnv.error_summary_group_map.clear()
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_FrAbs_N")
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_FrAbs_pi");
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_FrCEx_N");
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_FrCEx_pi");
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_FrElas_N");
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_FrElas_pi");
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_FrInel_N");
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_FrInel_pi");
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_FrPiProd_N");
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_FrPiProd_pi");
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_MFP_N");
mnv.error_summary_group_map["Genie_FSI"].push_back("genie_MFP_pi");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_AGKYxF1pi");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_AhtBY");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_BhtBY");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_CCQEPauliSupViaKF");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_CV1uBY");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_CV2uBY");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_EtaNCEL");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_MaCCQE");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_MaCCQEshape");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_MaNCEL");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_MaRES");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_MvRES");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_NormCCQE");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_NormCCRES");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_NormDISCC");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_NormNCRES");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_RDecBR1gamma");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_Rvn1pi");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_Rvn2pi");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_Rvn3pi");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_Rvp1pi");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_Rvp2pi");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_Theta_Delta2Npi");
mnv.error_summary_group_map["Genie_InteractionModel"].push_back("genie_VecFFCCQEshape");
mnv.DrawErrorSummary(mc_new)
#raw_input()
canvas1.Print("try1nuE.png")
