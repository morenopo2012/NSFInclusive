import os
import sys
import ROOT
from PlotUtils import MnvH1D,MnvH2D
_file = ROOT.TFile("Hists_Migration_t1_z26_Nu_v1_MnvG.root","READ")
hist = _file.Get("selected_mc_reco2d_Emu_Ehad")
for i in range(1,hist.GetNbinsX()+1):
    for j in range(1,hist.GetNbinsY()+1):
        for k in range(1,hist.GetNbinsZ()+1):
            cont = hist.GetBinContent(i,j,k)
            if cont==0:
                glob_bin = hist.FindBin(i,j,k)
                print glob_bin
