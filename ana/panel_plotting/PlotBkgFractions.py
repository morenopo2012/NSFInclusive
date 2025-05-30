import os
import sys
import PlotUtils


myplotter = PlotUtils.MnvPlotter()
myplotter.SetROOT6Palette(87)


f1 = ROOT.TFile("ScaledBkgCat_MichelSideBand.root")
f2 = ROOT.TFile("ScaledBkgCat_MichelSideBand.root")
f3 = ROOT.TFile("ScaledBkgCat_MichelSideBand.root")
f4 = ROOT.TFile("ScaledBkgCat_MichelSideBand.root")

