import ROOT
import PlotUtils
import os,sys


numfile = ROOT.TFile(sys.argv[1])
denfile = ROOT.TFile(sys.argv[2])

histname = "h_protonTndiff"

h_num = numfile.Get(histname)
h_den = denfile.Get(histname)

h_num_px = h_num.ProjectionX("allnum",0,-1)
h_den_px = h_den.ProjectionX("allden",0,-1)

h_num_px.Divide(h_num_px,h_den_px)

h_num_px.Draw()

raw_input("done")
