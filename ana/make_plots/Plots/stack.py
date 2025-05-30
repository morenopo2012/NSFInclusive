import ROOT
import os,sys
from PlotUtils import MnvH1D,MnvH2D,MnvVertErrorBand,MnvLatErrorBand,MnvPlotter
#the idea is to do the closure test.....
os.getenv("PERSISTENT")
#_file = ROOT.TFile("$HISTS/v1/Hists_Energy_MC_t1_z26_Nu_v1_.root","READ")


Emusig_list = ["Emu_trans","Emu_lowQ2Trans","Emu_lowQ2","Emu_lowW","Emu_dis","Emu_mc"]
#Emusig_list = ["Emu_mc"]
#sig_list = ["Emu_data"]
sig_list = ["Enu_trans","Enu_lowQ2Trans","Enu_lowQ2","Enu_lowW","Enu_dis","Enu_mc"]
#This script assumes it contains only the histograms that will go in the stack
#hist_list = ROOT.TObjArray()
#prefix_hist = "Enu_"
#signal_hist = prefix_hist+"Enu_"
global colors
colors = {1:ROOT.kBlack-6,
2:ROOT.kMagenta-6,
3:ROOT.kRed-6,
4:ROOT.kYellow-6,
5:ROOT.kWhite,
6:ROOT.kWhite,
7:ROOT.kWhite,
8:ROOT.kGreen-6,
9:ROOT.kTeal-6,
11:ROOT.kBlue-6,
12:ROOT.kMagenta-10,
13:ROOT.kRed-10,
14:ROOT.kYellow-10,
15:ROOT.kGray,
16:ROOT.kBlack,
17:ROOT.kBlack,
18:ROOT.kGreen-6,
19:ROOT.kTeal-6
}

def MakeHistPretty(histogram):
    thetitle = histogram.GetTitle()
    tot_entry = histogram.Integral()
    histogram.GetXaxis().SetTitle("Neutrino Energy (GeV)")
    histogram.GetYaxis().SetTitle("Entries/bin")
    histogram.GetXaxis().SetTitleSize(0.05)
    histogram.SetMarkerStyle(1)
    histogram.SetMarkerStyle(20)

    #histogram.SetMaximum(histogram.GetMaximum())
    #histogram.SetMinimum(1)

    histogram.GetYaxis().SetNdivisions(505)
    histogram.GetYaxis().SetTitleOffset(1.35)
    histogram.GetXaxis().SetTitleOffset(1.2)
    histogram.GetYaxis().CenterTitle()
    histogram.GetXaxis().CenterTitle()
    histogram.GetXaxis().SetNdivisions(505)
    histogram.Draw("hist")

def GetListOfHistograms(rfile):
    keys = rfile.GetListOfKeys()
    for key in keys:
        temp_name = key.GetName()
        if "mc" in temp_name:
            new_name = temp_name.replace("_mc","")
            print new_name
        else:
            continue

def EmptyErrorBands(vertnames,latnames,hist):
    for names in vertnames:
        if hist.HasVertErrorBand(names):
            hist.PopVertErrorBand(names)
    for names in latnames:
        if hist.HasLatErrorBand(names):
            hist.PopLatErrorBand(names)


def InitializeLegend():
    legend = ROOT.TLegend(0.55,0.5,0.9,0.87)
    legend.SetNColumns(2)
    legend.SetBorderSize(0)
    return legend

def EmptyErrorBands(vertnames,latnames,hist):
    for names in vertnames:
        if hist.HasVertErrorBand(names):
            hist.PopVertErrorBand(names)
    for names in latnames:
        if hist.HasLatErrorBand(names):
            hist.PopLatErrorBand(names)

leg = InitializeLegend()


 
it_list = [2,3,11,18,4,1]
mc = ROOT.THStack("newstack","")
rfile = ROOT.TFile("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_.root","READ")


canvas1 = ROOT.TCanvas("c","c",200,10,600,480)
#canvas = ROOT.TCanvas()
canvas1.cd()
for i in range(len(sig_list)):
    temp_hist = rfile.Get(sig_list[i])
#for key in rfile.GetListOfKeys(len(sig_list)):
    #ss = i.GetName()
    #temp_hist = _file.Get(key.GetName())
    temp_hist.SetDirectory(0)
    temp_hist.SetFillColor(colors[it_list[i]])
    ss=sig_list[i]
    leg.AddEntry(temp_hist,ss,"f")
    temp_hist.SetTitle(sig_list[i])
    #MakeHistPretty(temp_hist)
    mc.Add(temp_hist)
    mc.Draw("hist")

#if you have all kinds of histograms in your rootfile you need to do it one by one
mc.GetXaxis().SetTitleSize(0.05)
mc.GetYaxis().SetTitleOffset(1.35)
mc.GetXaxis().SetTitleOffset(1.2)
mc.GetYaxis().CenterTitle()
mc.GetXaxis().CenterTitle()
#mc.SetMarkerStyle(1)
#mc.SetMarkerStyle(20)
mc.GetXaxis().SetTitle("Neutrino Energy (GeV)")
mc.GetYaxis().SetTitle("Entries/bin")
leg.Draw()
raw_input()
canvas1.Print("PlotStack.png") 


leg1 = InitializeLegend()
mc1 = ROOT.THStack("newstack","")
canvas2 = ROOT.TCanvas("c","c",200,10,600,480)
#canvas = ROOT.TCanvas()
canvas2.cd()
for i in range(len(Emusig_list)):
    temp_hist = rfile.Get(Emusig_list[i])
#for key in rfile.GetListOfKeys(len(sig_list)):
    #ss = i.GetName()
    #temp_hist = _file.Get(key.GetName())
    temp_hist.SetDirectory(0)
    temp_hist.SetFillColor(colors[it_list[i]])
    ss=Emusig_list[i]
    leg1.AddEntry(temp_hist,ss,"f")
    temp_hist.SetTitle(Emusig_list[i])
    #MakeHistPretty(temp_hist)
    mc1.Add(temp_hist)
    mc1.Draw("hist")

#if you have all kinds of histograms in your rootfile you need to do it one by one
mc1.GetXaxis().SetTitleSize(0.05)
mc1.GetYaxis().SetTitleOffset(1.35)
mc1.GetXaxis().SetTitleOffset(1.2)
mc1.GetYaxis().CenterTitle()
mc1.GetXaxis().CenterTitle()
#mc.SetMarkerStyle(1)
#mc.SetMarkerStyle(20)
mc1.GetXaxis().SetTitle("Muon Energy (GeV)")
mc1.GetYaxis().SetTitle("Entries/bin")
leg1.Draw()
raw_input()
canvas2.Print("EmuPlotStack.png") 

