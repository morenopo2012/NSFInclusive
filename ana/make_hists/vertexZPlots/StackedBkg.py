import ROOT
import os,sys
from ROOT import *
from ROOT import PlotUtils
from ROOT import TLegend
from ROOT import THStack
from ROOT import THStack
infile= ROOT.TFile("Hists_EventSelection_t99_z99_Nu_minervame1L.root")
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()


mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()
mcScale =  dataPOT/mcPOT



water = infile.Get("selected_mc_reco_water_vtxz")
carbon = infile.Get("selected_mc_reco_carbon_vtxz")
iron = infile.Get("selected_mc_reco_iron_vtxz")
lead = infile.Get("selected_mc_reco_lead_vtxz")
scint = infile.Get("selected_mc_reco_scintillator_vtxz")
data_hist = infile.Get("selected_data_reco_vtxz")

# do not need to be width scale, all bins are the same
#water.Scale(mcScale)

carbon.Scale(mcScale)
iron.Scale(mcScale)
lead.Scale(mcScale)
scint.Scale(mcScale)

mcColors = MnvColors.GetColors(MnvColors.kOkabeItoPalette)
water.SetLineColor(mcColors[4])
carbon.SetLineColor(mcColors[5])
iron.SetLineColor(mcColors[2])
lead.SetLineColor(12)
    #mcColors[7])
scint.SetLineColor(mcColors[3])
water.SetFillColor(mcColors[4])
carbon.SetFillColor(mcColors[5])
iron.SetFillColor(mcColors[2])
lead.SetFillColor(12)
    #mcColors[7])
scint.SetFillColor(mcColors[3])
water.Scale(water.GetNormBinWidth(),"width")
carbon.Scale(carbon.GetNormBinWidth(),"width")
lead.Scale(lead.GetNormBinWidth(),"width")
iron.Scale(iron.GetNormBinWidth(),"width")
scint.Scale(scint.GetNormBinWidth(),"width")
data_hist.Scale(data_hist.GetNormBinWidth(),"width")
#water.Scale(1/water.Integral())
#carbon.Scale(1/carbon.Integral())
#lead.Scale(1/lead.Integral())
#iron.Scale(1/iron.Integral())
#scint.Scale(1/scint.Integral())
#data_hist.Scale(1/data_hist.Integral())


data_hist.SetMarkerStyle(20)
data_hist.SetMarkerSize(1)
data_hist.SetMarkerColor(1)
data_hist.SetLineWidth(1)
data_hist.SetLineStyle(1)
data_hist.SetLineColor(1)
data_hist.SetFillColor(0)


#scint.Draw("HIST")
lead.Draw("HIST")
iron.Draw("HIST SAME")
carbon.Draw("HIST SAME")
water.Draw("HIST SAME")


data_hist.Draw("HIST SAME p E1 X0") # for error bars, suppressed error bars along X
stack = THStack("stack","stack")
stack.Add(scint)
stack.Add(carbon)
stack.Add(water)
stack.Add(iron)
stack.Add(lead)
#stack.Add(data_hist)
stack.Draw("HIST")


data_hist.Draw("HIST SAME p E1 X0") # for error bars, suppressed error bars along X
stack.GetXaxis().SetTitle("Reconstructed vertex Z (cm)")
stack.GetYaxis().SetTitle("Events/cm")
stack.GetXaxis().CenterTitle()
stack.GetYaxis().CenterTitle()
stack.GetYaxis().SetTitleOffset(1.1)
stack.GetXaxis().SetRangeUser(440, 590)
stack.SetMaximum(stack.GetMaximum()*1.44)
mnv.AddHistoTitle("ML vertexing", 0.05, 1)
#mnv.WritePreliminary(0.38, 0.76, 0.035, True)
#mnv.WriteNorm("POT-Normalized", 0.29, 0.81, 0.035 )
mnv.WritePreliminary(0.38, 0.73, 0.035, True)
mnv.WriteNorm("POT-Normalized", 0.29, 0.77, 0.035 )
#mnv.WritePreliminary(0.36, 0.82, 0.035, True)
#legend = TLegend(0.60,0.64,0.90,0.89)
#legend = TLegend(0.17,0.85,0.88,0.89)
legend = TLegend(0.20,0.80,0.75,0.89)
legend.SetNColumns(3)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextSize(0.035)
#legend.AddEntry(data_hist, " Data", "ep")
legend.AddEntry(water, " Water", "fl")
legend.AddEntry(carbon, " Carbon", "fl")
legend.AddEntry(iron, " Iron", "fl")
legend.AddEntry(lead, " Lead", "fl")
legend.AddEntry(scint, " Scintillator", "fl")
legend.SetTextFont(42)
legend.Draw()
#canvas1.SetLogy()
canvas1.Modified()
canvas1.Print("ReconstructedVertexZ_TBV_data_POTNorm_water_NEW_NOLOG.png")
raw_input("Done")
