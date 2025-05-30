import ROOT
import os,sys,math
from PlotUtils import *
import array
_file = ROOT.TFile(str(sys.argv[1]),"READ")
mnv = MnvPlotter()
mnv.ApplyStyle(8)
#_migration = _file.Get("selected_mc_response1d_Emu_migration")
#_migration = _file.Get("selected_mc_response2d_Emu_Ehad_migration")
#_migration = _file.Get("selected_mc_response1d_Ehad_migration")

#_migration = _file.Get("selected_mc_response2d_Emu_Ehad_migration")
#_migration = _file.Get("selected_mc_response2d_W_Q2_migration")
#_migration = _file.Get("selected_mc_response2d_x_Q2_migration")
_migration = _file.Get("selected_mc_response2d_x_y_migration")

#_migration = _file.Get("h_enu_migration_qelike")
#_migration = _file.Get("h_pzmu_ptmu_qelike_migration")
#_migration = _file.Get("h_ptmu_qelike_migration")
_normmigration = _migration.Clone()
#now the row normalized migration....
nbinsy = _migration.GetNbinsY()
nbinsx = _migration.GetNbinsX()
for i in range(0,nbinsy+1):
    row_norm = 0.0
    for j in range(0,nbinsx+1):
        row_norm += _migration.GetBinContent(j,i)
    for j in range(0,nbinsy+1):
        _cont = _normmigration.GetBinContent(j,i)
        if row_norm!=0.0:
            _normmigration.SetBinContent(j,i,_cont/row_norm)
mnv = MnvPlotter()
mnv.SetCorrelationPalette()
canvas = ROOT.TCanvas("c","c",750,750)
canvas.cd()
_normmigration.SetMaximum(1.0)
#_normmigration.GetXaxis().SetTitle("Reconstructed Muon P_{||} per Muon P_{t} bins")
#_normmigration.GetYaxis().SetTitle("True Muon P_{||} per Muon P_{t} bins")
### Along X axis --> Reconstructed X per Y variable
## Along Y axis --> True X per Y variable

#_normmigration.GetXaxis().SetTitle("Reco E_{#mu} per Ehad bins")
#_normmigration.GetYaxis().SetTitle("True E_{#mu} per Ehad bins")
#_normmigration.GetXaxis().SetTitle("Reco W per Q2 Bins")
#_normmigration.GetYaxis().SetTitle("True W per Q2 Bins")
#_normmigration.GetXaxis().SetTitle("Reco x per Q2 Bins")
#_normmigration.GetYaxis().SetTitle("True x per Q2 Bins")
_normmigration.GetXaxis().SetTitle("Reco x per y Bins")
_normmigration.GetYaxis().SetTitle("True x per y Bins")

#_normmigration.GetYaxis().SetTitle("Reconstructed E_{#nu}(GeV)")
#_normmigration.GetXaxis().SetTitle("True E_{#nu}(GeV)")
_normmigration.Draw("colz")
raw_input()
