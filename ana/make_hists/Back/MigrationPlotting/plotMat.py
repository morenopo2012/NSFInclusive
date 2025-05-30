import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle

def colNormalize(hist):
        norm_hist = hist.Clone()
        norm_hist.Reset()

        nXBins = hist.GetNbinsX()+2 # to add overflow + underflow
        nYBins = hist.GetNbinsY()+2 # to add overflow + underflow
        
        for i in range(0, nXBins): # column
                col_integral = 0
                for j in range(0,nYBins): # row
                        col_integral += hist.GetBinContent(i,j)
                for j in range(0,nYBins):
                        if col_integral==0: norm_hist.SetBinContent(i,j,0)
                        else: norm_hist.SetBinContent(i,j, hist.GetBinContent(i,j)*100/col_integral)
        return norm_hist


def rowNormalize(hist):
        norm_hist = hist.Clone()
        norm_hist.Reset()

        nXBins = hist.GetNbinsX()+2 # to add overflow + underflow
        nYBins = hist.GetNbinsY()+2 # to add overflow + underflow
        
        for i in range(0, nYBins): # column
                row_integral = 0
                for j in range(0,nXBins): # row
                        row_integral += hist.GetBinContent(j,i)
                for j in range(0,nXBins):
                        if row_integral==0: norm_hist.SetBinContent(j,i,0)
                        else: norm_hist.SetBinContent(j,i, hist.GetBinContent(j,i)*100/row_integral)
        return norm_hist


targetID = sys.argv[1] 
targetZ = sys.argv[2]
playlist = sys.argv[3]
infile= ROOT.TFile("/minerva/app/users/zdar/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive/ana/make_hists/Hists_Migration_t%s_z%s_AntiNu_%s.root"%(targetID,targetZ,playlist),"READ")
mnv = PlotUtils.MnvPlotter()
canvas1 = ROOT.TCanvas() 
ROOT.TH1.AddDirectory(False)

vars = ["Emu", "Ehad"]  
for myvariable in vars:

 mc_hist = infile.Get("response1d_%s_migration"%myvariable)

# to add overflow + underflow
 mc_hist.GetXaxis().SetRangeUser(-1,mc_hist.GetNbinsX()+1)
 mc_hist.GetYaxis().SetRangeUser(-1,mc_hist.GetNbinsY()+1)

 gStyle. SetPalette(1)
 gStyle.SetOptTitle(0)
 gStyle.SetPaintTextFormat("2.0f")

 if myvariable == "Enu":
        mc_hist.GetXaxis().SetTitle("Reconstructed Antineutrino Energy Bins")
        mc_hist.GetYaxis().SetTitle("True Antineutrino Energy Bins")
 if myvariable == "Ehad":
        mc_hist.GetXaxis().SetTitle("Reconstructed Hadronic Energy Bins")
        mc_hist.GetYaxis().SetTitle("True Hadronic Energy Bins")
 if myvariable == "x":
        mc_hist.GetXaxis().SetTitle("Reconstructed Bjorken x Bins")
        mc_hist.GetYaxis().SetTitle("True Bjorken x Bins")
 if myvariable == "Emu":
        mc_hist.GetXaxis().SetTitle("Reconstructed Muon Energy Bins")
        mc_hist.GetYaxis().SetTitle("True Muon Energy Bins")

 mc_hist.GetXaxis().CenterTitle()
 mc_hist.GetXaxis().SetTitleFont(42)
 mc_hist.GetXaxis().SetTitleSize(0.05)

 mc_hist.GetYaxis().SetTitleOffset(1)
 mc_hist.GetYaxis().CenterTitle()
 mc_hist.GetYaxis().SetTitleFont(42)
 mc_hist.GetYaxis().SetTitleSize(0.05)
 
 mc_hist.GetZaxis().CenterTitle()
 mc_hist.GetZaxis().SetRangeUser(0,1000)
 mc_hist.GetZaxis().SetTitleFont(42)
 mc_hist.GetZaxis().SetTitleOffset(1.1)
 mc_hist.GetZaxis().SetTitleSize(0.045)
 mc_hist.SetMarkerSize(1)
 
 if targetZ == "26" and targetID == "1":
  mnv.AddHistoTitle("Iron of target 1", 0.042, 1)
 if targetZ == "26" and targetID == "2":
  mnv.AddHistoTitle("Iron of target 2", 0.042, 1)
 if targetZ == "26" and targetID == "3":
  mnv.AddHistoTitle("Iron of target 3", 0.042, 1)
 if targetZ == "26" and targetID == "5":
  mnv.AddHistoTitle("Iron of target 5", 0.042, 1)
 if targetZ == "06" and targetID == "3":
  mnv.AddHistoTitle("Carbon of target 3", 0.042, 1)
 if targetZ == "82" and targetID == "1":
  mnv.AddHistoTitle("Lead of target 1", 0.042, 1)
 if targetZ == "82" and targetID == "2":
  mnv.AddHistoTitle("Lead of target 2", 0.042, 1)
 if targetZ == "82" and targetID == "3":
  mnv.AddHistoTitle("Lead of target 3", 0.042, 1)
 if targetZ == "82" and targetID == "4":
  mnv.AddHistoTitle("Lead of target 4", 0.042, 1)
 if targetZ == "82" and targetID == "5":
  mnv.AddHistoTitle("Lead of target 5", 0.042, 1)
 if targetZ == "82" and targetID == "14":
  mnv.AddHistoTitle("Tracker modules 27-32", 0.042, 1)
 if targetZ == "82" and targetID == "24":
  mnv.AddHistoTitle("Tracker modules 33-38", 0.042, 1)
 if targetZ == "82" and targetID == "34":
  mnv.AddHistoTitle("Tracker modules 39-44", 0.042, 1)
 if targetZ == "82" and targetID == "44":
  mnv.AddHistoTitle("Tracker modules 45-50", 0.042, 1)
 if targetZ == "82" and targetID == "54":
  mnv.AddHistoTitle("Tracker modules 51-56", 0.042, 1)
 if targetZ == "82" and targetID == "64":
  mnv.AddHistoTitle("Tracker modules 57-62", 0.042, 1)
 if targetZ == "82" and targetID == "74":
  mnv.AddHistoTitle("Tracker modules 63-68", 0.042, 1)
 if targetZ == "82" and targetID == "84":
  mnv.AddHistoTitle("Tracker modules 69-74", 0.042, 1)
 if targetZ == "82" and targetID == "94":
  mnv.AddHistoTitle("Tracker modules 75-80", 0.042, 1)
 # SIMPLE OCCUPANCY
 mc_hist.Draw("COLZTEXT")
 mnv.SetROOT6Palette(109)
 mc_hist.GetZaxis().SetTitle("Number of Events")
 canvas1.Print("Mig_t%s_z%s_%s_Occupancy.png"%(targetID,targetZ,myvariable))
 
 
 # COLUMN NORMALIZED
 mc_hist.GetZaxis().SetRangeUser(0,100)
 mc_hist.GetZaxis().SetTitle("Column Normalized Event Rate (%)")
 colnorm_hist = colNormalize(mc_hist)
 colnorm_hist.Draw("COLZTEXT")
 canvas1.Print("Mig_t%s_z%s_%s_ColumnNorm.png"%(targetID,targetZ,myvariable))

 # ROW NORMALIZED
 mc_hist.GetZaxis().SetTitle("Row Normalized Event Rate (%)")
 rownorm_hist = rowNormalize(mc_hist)
 rownorm_hist.Draw("COLZTEXT")
 canvas1.Print("Mig_t%s_z%s_%s_RowNorm.png"%(targetID,targetZ,myvariable)) 
raw_input("Done")
