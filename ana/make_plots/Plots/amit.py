import ROOT
import os,sys,math
from PlotUtils import *
import array
#_file = ROOT.TFile("/minerva/app/users/zdar/cmtuser/Minerva_v21r1p1_MasterAnaDev/Ana/NSFNukeCCInclusive/ana/make_hists/Hists_EventSelection_t14_z82_Nu_v1_.root","READ")
_file = ROOT.TFile("/minerva/app/users/zdar/cmtuser/Minerva_v21r1p1_MasterAnaDev/Ana/NSFNukeCCInclusive/ana/make_hists/Hists_EventSelection_t14_z82_Nu_v1_.root","READ")
#_file = ROOT.TFile("/minerva/data/users/zdar/NukeHists/MuonKludged/minervame1L/Hists_Energy_MC_t14_z82_Nu_MuonKludged.root","READ")
mc = _file.Get("selected_mc_reco_Emu")
#mc = _file.Get("sample_MC_Emu_t14_z82")
#mc = _file.Get("Emu_mc")
vertnames = mc.GetVertErrorBandNames()
for name in vertnames:
    #mc_x = mc.ProjectionX()
    #dat_x = dat_bkgsub.ProjectionX()
    _mc = mc.GetVertErrorBand(name).GetErrorBand(True)
    _mc.SetTitle(name+" Before Bkg sub.")
    _mc.SetLineColor(2)
    canvas = ROOT.TCanvas()
    canvas.cd()
    _mc.SetMinimum(0.0)
    _mc.Draw("hist")
    #print _dat.GetMaximum(),_mc.GetMaximum()
    _name = name.replace(" ","")
    canvas.Print(_name+".png")
