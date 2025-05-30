import ROOT
from ROOT import PlotUtils
import numpy as np

# Load the ROOT file
f = ROOT.TFile.Open("Hists_EventSelection_t95_z95_Nu_minervame1G.root")

# Define histograms
mc_total = f.Get("m_selected_mc_reco_QE_Eavailable_q3").Clone("mc_total")
mc_QE    = f.Get("m_selected_mc_reco_QE_Eavailable_q3")
mc_RES   = f.Get("m_selected_mc_reco_RES_Eavailable_q3")
mc_DIS   = f.Get("m_selected_mc_reco_DIS_Eavailable_q3")
mc_2p2h  = f.Get("m_selected_mc_reco_2p2h_Eavailable_q3")
mc_Other = f.Get("m_selected_mc_reco_OtherIT_Eavailable_q3")
data     = f.Get("h_data_Eavailable_q3")

# POT Normalization
data_pot = f.Get("DataPOT").GetVal()
mc_pot = f.Get("MCPOT").GetVal()
pot_ratio = data_pot / mc_pot

# Rebinning info
q3_bins = [0.0, 0.2, 0.3, 0.4, 0.6, 0.9, 1.2]  # GeV

# Canvas for plotting
canvas = ROOT.TCanvas("canvas", "", 1200, 1800)
canvas.Divide(2, 3)

for i in range(len(q3_bins) - 1):
    canvas.cd(i+1)

    q3_min, q3_max = q3_bins[i], q3_bins[i+1]

    # Get bin ranges for q3 axis
    q3_axis = data.GetYaxis()
    bin_min = q3_axis.FindBin(q3_min + 1e-4)
    bin_max = q3_axis.FindBin(q3_max - 1e-4)

    # Project each component onto X-axis (Eavail) for this q3 bin
    proj_data = data.ProjectionX(f"proj_data_{i}", bin_min, bin_max)
    proj_QE = mc_QE.ProjectionX(f"proj_QE_{i}", bin_min, bin_max)
    proj_RES = mc_RES.ProjectionX(f"proj_RES_{i}", bin_min, bin_max)
    proj_DIS = mc_DIS.ProjectionX(f"proj_DIS_{i}", bin_min, bin_max)
    proj_2p2h = mc_2p2h.ProjectionX(f"proj_2p2h_{i}", bin_min, bin_max)
    proj_Other = mc_Other.ProjectionX(f"proj_Other_{i}", bin_min, bin_max)

    # Scale MC by POT
    for h in [proj_QE, proj_RES, proj_DIS, proj_2p2h, proj_Other]:
        h.Scale(pot_ratio)

    # Stack MC histograms
    stack = ROOT.THStack(f"stack_{i}", "")
    proj_QE.SetLineColor(ROOT.kBlue)
    proj_RES.SetLineColor(ROOT.kViolet)
    proj_DIS.SetLineColor(ROOT.kGreen+2)
    proj_2p2h.SetLineColor(ROOT.kOrange+7)
    proj_Other.SetLineColor(ROOT.kGray+1)

    for h in [proj_Other, proj_2p2h, proj_DIS, proj_RES, proj_QE]:
        h.SetFillColor(h.GetLineColor())
        stack.Add(h)

    # Total MC
    totalMC = proj_QE.Clone(f"totalMC_{i}")
    for h in [proj_RES, proj_DIS, proj_2p2h, proj_Other]:
        totalMC.Add(h)

    totalMC.SetLineColor(ROOT.kRed+1)
    totalMC.SetLineWidth(2)

    # Draw stack and total
    proj_data.SetMarkerStyle(20)
    proj_data.SetMarkerSize(0.8)
    proj_data.SetLineColor(ROOT.kBlack)
    proj_data.SetTitle(f"{q3_min:.2f} < q₃ < {q3_max:.2f} GeV")

    stack.Draw("HIST")
    totalMC.Draw("HIST SAME")
    proj_data.Draw("E1 SAME")

    # Make ratio plot
    ratio = proj_data.Clone(f"ratio_{i}")
    ratio.Divide(totalMC)
    ratio.SetTitle("")
    ratio.GetYaxis().SetTitle("Data / MC")
    ratio.GetYaxis().SetRangeUser(0, 2)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(0.6)
    ratio.SetLineColor(ROOT.kBlack)

    pad = ROOT.TPad(f"pad_ratio_{i}", "", 0, 0, 1, 0.3)
    canvas.cd(i+1)
    pad.SetTopMargin(0.01)
    pad.SetBottomMargin(0.3)
    pad.Draw()
    pad.cd()
    ratio.Draw("E1")

canvas.SaveAs("Eavailable_vs_q3_panels.pdf")

