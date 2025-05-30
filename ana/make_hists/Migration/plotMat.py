import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle
#infile= ROOT.TFile("")i
#infile = sys.argv[1]
canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
mnv = PlotUtils.MnvPlotter()

gStyle. SetPalette(1)
gStyle.SetOptTitle(0)


# Open the ROOT file
#infile = ROOT.TFile.Open("Migration_minervame1L_t99_z99_nosys.root")
infile = ROOT.TFile.Open("Migration_minervame1G_t95_z95_nosys.root") #Marvin's thesis
#myvariable = sys.argv[2]

myvariable = "Eavailable_q3"

#mc_hist = infile.Get("response2d_%s_migration"%myvariable)
mc_hist = infile.Get("selected_mc_response2d_%s_migration"%myvariable)
#print(mc_hist.GetNbinsX())
#print(mc_hist.GetNbinsY())
#mc_hist.Draw("COLZ TEXTnn")
mnv.ApplyStyle(6)
mnv.DrawNormalizedMigrationHistogram(mc_hist, True, False, True, True)
'''
void MnvPlotter::DrawNormalizedMigrationHistogram(
        const TH2D* h_migration,
        const bool drawAsMatrix,
        const bool coarseContours, /* = false */
        const bool includeFlows, /* = true */
        const bool noText /* = false */
        )
{
'''
# ROW NORMALIZED
#mnv.AddHistoTitle("Fe Target 2", 0.035, 1)  
canvas1.Modified()
mnv.SetROOT6Palette(113)

keys = canvas1.GetListOfPrimitives();
for k in keys:
    if(k.ClassName().find("TH2D")!=-1): # to change axis titles
        if myvariable == "Emu_Ehad":
                k.GetXaxis().SetTitle("Reconstructed Emu  per Ehad bins")
                k.GetYaxis().SetTitle("True Emu per Ehad bins")
        if myvariable == "Enu_Ehad":
                k.GetXaxis().SetTitle("Reconstructed Enu per Ehad Bins")
                k.GetYaxis().SetTitle("True Enu per Ehad Bins")
        if myvariable == "x_Q2":
                k.GetXaxis().SetTitle("Reconstructed x_{Bj} per Q^{2} Bins")
                k.GetYaxis().SetTitle("True x_{Bj} per Q^{2} Bins")
        if myvariable == "W_Q2":
                k.GetXaxis().SetTitle("Reconstructed W per Q^{2} Bins")
        if myvariable == "Eavailable_q3":
                k.GetXaxis().SetTitle("Reconstructed E_{av} per q_{3} Bins")
                k.GetYaxis().SetTitle("True E_{av} per q_{3} Bins")

   #     if myvariable == "xfine":
    #            k.GetXaxis().SetTitle("Reconstructed x_Bj Bins")
     #          k.GetYaxis().SetTitle("True x_Bj Bins")

        k.GetXaxis().CenterTitle()
        k.GetXaxis().SetTitleFont(42)
        k.GetXaxis().SetTitleSize(0.05)

        k.GetYaxis().SetTitleOffset(1)
        k.GetYaxis().CenterTitle()
        k.GetYaxis().SetTitleFont(42)
        k.GetYaxis().SetTitleSize(0.05)

        k.GetZaxis().SetTitle("Row Normalized Event Rate (%)")
        k.GetZaxis().CenterTitle()
        k.GetZaxis().SetRangeUser(0,100)
        k.GetZaxis().SetTitleFont(42)
        k.GetZaxis().SetTitleOffset(1.1)
        k.GetZaxis().SetTitleSize(0.045)
    
canvas1.SetLeftMargin(0.12)
canvas1.SetBottomMargin(0.12)
canvas1.SetRightMargin(0.16)
canvas1.Print("MigrationMatrix_%s.png"%myvariable)

raw_input("Done")
