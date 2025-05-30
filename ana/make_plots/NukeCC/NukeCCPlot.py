import ROOT
import os,sys
from PlotUtils import MnvH1D,MnvH2D,MnvPlotter
#_fnsf = ROOT.TFile("/minerva/data/users/afilkins/NukeHists/MuonKludged/minervame1L/Hists_Energy_MC_t14_z82_Nu_MuonKludged.root","READ")
_fnsf = ROOT.TFile("/pnfs/minerva/persistent/users/jleclerc/NukeHists/MuonKludged/minervame1L/Hists_Energy_MC_t14_z82_Nu_MuonKludged.root","READ")
#_fnsf = ROOT.TFile("Hists_Energy_v1_t14_z82_Nu_v1_.root","READ")
#_fnsf = ROOT.TFile("Hists_EventSelection_t1_z26_Nu_v1_.root","READ")
#_fnsf = ROOT.TFile("Hists_EventSelection_t14_z82_Nu_v1_.root","READ")
#mc_name = "Enu_mc"
#mc_name = "selected_mc_reco_Emu"
mc_name = "sample_MC_Emu_t14_z82"
mc_new = _fnsf.Get(mc_name)


mnv = MnvPlotter()
canvas1 = ROOT.TCanvas()
mnv.error_summary_group_map.clear()
#mnv.error_summary_group_map["Muon"].push_back("Muon_Energy_Resolution");
#mnv.error_summary_group_map["Muon"].push_back("Muon_Energy_MINOS");
#mnv.error_summary_group_map["MINOS"].push_back("MINOS_Reconstruction_Efficiency");
mnv.error_summary_group_map["Det Res"].push_back("Vertex Rec.");
mnv.error_summary_group_map["Det Res"].push_back("Muon Energy Rec.");
mnv.error_summary_group_map["Det Res"].push_back("Hadronic Energy Rec.");
mnv.error_summary_group_map["Det Res"].push_back("Birks' Parameter");
mnv.error_summary_group_map["Flux"].push_back("Flux");
mnv.error_summary_group_map["MASS"].push_back("CH Mass");
mnv.error_summary_group_map["MASS"].push_back("CMass");
mnv.error_summary_group_map["MASS"].push_back("FeMass");
mnv.error_summary_group_map["MASS"].push_back("PbMass");
mnv.error_summary_group_map["MASS"].push_back("Normalization");
mnv.error_summary_group_map["GENIE"].push_back("Non Resonant Pion");
mnv.error_summary_group_map["RPA"].push_back("RPA-Model");
mnv.error_summary_group_map["GENIE"].push_back("2p2h-Model");
mnv.error_summary_group_map["GENIE"].push_back("AGKY Model");
mnv.error_summary_group_map["GENIE"].push_back("FrAbs_N");
mnv.error_summary_group_map["GENIE"].push_back("FrAbs_pi");
mnv.error_summary_group_map["GENIE"].push_back("FrCEx_N");
mnv.error_summary_group_map["GENIE"].push_back("FrCEx_pi");
mnv.error_summary_group_map["GENIE"].push_back("FrElas_N");
mnv.error_summary_group_map["GENIE"].push_back("FrElas_pi");
mnv.error_summary_group_map["GENIE"].push_back("FrInel_N");
mnv.error_summary_group_map["GENIE"].push_back("FrInel_pi");
mnv.error_summary_group_map["GENIE"].push_back("FrPiProd_N");
mnv.error_summary_group_map["GENIE"].push_back("FrPiProd_pi");
mnv.error_summary_group_map["GENIE"].push_back("MFP_N");
mnv.error_summary_group_map["GENIE"].push_back("MFP_pi");
mnv.error_summary_group_map["GENIE"].push_back("AGKYxF1pi");
mnv.error_summary_group_map["GENIE"].push_back("AhtBY");
mnv.error_summary_group_map["GENIE"].push_back("BhtBY");
mnv.error_summary_group_map["GENIE"].push_back("CCQEPauliSupViaKF");
mnv.error_summary_group_map["GENIE"].push_back("CV1uBY");
mnv.error_summary_group_map["GENIE"].push_back("CV2uBY");
mnv.error_summary_group_map["GENIE"].push_back("EtaNCEL");
mnv.error_summary_group_map["GENIE"].push_back("MaCCQE");
mnv.error_summary_group_map["GENIE"].push_back("NormCCQE");
mnv.error_summary_group_map["GENIE"].push_back("MaNCEL");
mnv.error_summary_group_map["GENIE"].push_back("MaRES");
mnv.error_summary_group_map["GENIE"].push_back("MvRES");
mnv.error_summary_group_map["GENIE"].push_back("NormDISCC");
mnv.error_summary_group_map["GENIE"].push_back("NormNCRES");
mnv.error_summary_group_map["GENIE"].push_back("RDecBR1gamma");
mnv.error_summary_group_map["GENIE"].push_back("Rvn1pi");
mnv.error_summary_group_map["GENIE"].push_back("Rvn2pi");
mnv.error_summary_group_map["GENIE"].push_back("Rvn3pi");
mnv.error_summary_group_map["GENIE"].push_back("Rvp1pi");
mnv.error_summary_group_map["GENIE"].push_back("Rvp2pi");
mnv.error_summary_group_map["GENIE"].push_back("Theta_Delta2Npi");
mnv.error_summary_group_map["GENIE"].push_back("VecFFCCQEshape");
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
