import ROOT
import os,sys
from PlotUtils import MnvH1D,MnvH2D,MnvPlotter
_fnsf = ROOT.TFile("/minerva/app/users/zdar/cmtuser/Minerva_v21r1p1_MasterAnaDev/Ana/NSFNukeCCInclusive/ana/make_hists/Hists_Migration_t1_z26_Nu_v1_MnvGenie.root","READ")
#_fnsf = ROOT.TFile("/minerva/data/NSF_Validation/referenceHists/NSF_NukeCC_mc_minervame1L_2020-04-13.root","READ")
#mc_name = "Emu_mc"
mc_name = "selected_mc_reco_Emu"
#mc_name = "selected_data_reco_Emu"
#mc_name = "h_inclusive_Emu"
mc_new = _fnsf.Get(mc_name)


mnv = MnvPlotter()
canvas1 = ROOT.TCanvas()
mnv.error_summary_group_map.clear()
mnv.error_summary_group_map["Muon"].push_back("Muon_Energy_Resolution");
mnv.error_summary_group_map["Muon"].push_back("Muon_Energy_MINOS");
mnv.error_summary_group_map["Muon"].push_back("Muon_Energy_MINERvA");
mnv.error_summary_group_map["MINOS"].push_back("MINOS_Reconstruction_Efficiency");
mnv.error_summary_group_map["Det Res"].push_back("BeamAngleX");
mnv.error_summary_group_map["Det Res"].push_back("BeamAngleY");
mnv.error_summary_group_map["Flux"].push_back("Flux");
mnv.error_summary_group_map["RPA"].push_back("RPA_HighQ2");
mnv.error_summary_group_map["RPA"].push_back("RPA_LowQ2");
mnv.error_summary_group_map["GENIE"].push_back("Low_Recoil_2p2h_Tune");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrAbs_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrAbs_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrCEx_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrCEx_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrElas_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrElas_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrInel_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrInel_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrPiProd_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_FrPiProd_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_MFP_N");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_MFP_pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_AGKYxF1pi");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_AhtBY");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_BhtBY");
mnv.error_summary_group_map["GENIE"].push_back("GENIE_CCQEPauliSupViaKF");
#mnv.error_summary_group_map["GENIE"].push_back("CV1uBY");
#mnv.error_summary_group_map["GENIE"].push_back("CV2uBY");
#mnv.error_summary_group_map["GENIE"].push_back("EtaNCEL");
#mnv.error_summary_group_map["GENIE"].push_back("MaCCQE");
#mnv.error_summary_group_map["GENIE"].push_back("NormCCQE");
#mnv.error_summary_group_map["GENIE"].push_back("MaCCQEshape");
#mnv.error_summary_group_map["GENIE"].push_back("MaNCEL");
#mnv.error_summary_group_map["GENIE"].push_back("MaRES");
#mnv.error_summary_group_map["GENIE"].push_back("MvRES");
#mnv.error_summary_group_map["GENIE"].push_back("NormDISCC");
#mnv.error_summary_group_map["GENIE"].push_back("NormCCRES");
#mnv.error_summary_group_map["GENIE"].push_back("NormNCRES");
#mnv.error_summary_group_map["GENIE"].push_back("RDecBR1gamma");
#mnv.error_summary_group_map["GENIE"].push_back("Rvn1pi");
#mnv.error_summary_group_map["GENIE"].push_back("Rvn2pi");
#mnv.error_summary_group_map["GENIE"].push_back("Rvn3pi");
#mnv.error_summary_group_map["GENIE"].push_back("Rvp1pi");
#mnv.error_summary_group_map["GENIE"].push_back("Rvp2pi");
#mnv.error_summary_group_map["GENIE"].push_back("Theta_Delta2Npi");
#mnv.error_summary_group_map["GENIE"].push_back("VecFFCCQEshape");
#mnv.DrawErrorSummary(mc_new, "TR", 0, 1, 0.0, 0, "MinosEfficiency");
mnv.DrawErrorSummary(mc_new);
#mnv.DrawDataMCWithErrorBand(mc_new);
#mnv.DrawErrorSummary(mc_new, "TR", "false", "true", 0.0)
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "GENIE")
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "GENIE")
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "Flux")
#mnv.DrawErrorSummary(mc_new, "TR", 1, 1, 0.0, 0, "VecFFCCQEshape")
canvas1.Print("abc.png")
input("Press enter to end program")
