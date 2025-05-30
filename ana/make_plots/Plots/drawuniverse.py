import ROOT
from PlotUtils import MnvH1D, MnvLatErrorBand,MnvVertErrorBand
import os,sys
#open the file
if len(sys.argv)<3:
    print "checksystematics <filename> <histname>"
    sys.exit(1)

filename  = str(sys.argv[1])
hist_name = str(sys.argv[2])

rootfile = ROOT.TFile("Hists_EventSelection_t14_z82_Nu_v1_.root","READONLY")
histo = rootfile.Get("selected_mc_reco_Emu").Clone()  #beamfit_reweight_function

#latnames=histo.GetLatErrorBandNames()
vertnames=histo.GetVertErrorBandNames()
canvas= ROOT.TCanvas("canavas","canvas",750,750)
canvas.cd()
for names in vertnames:
    #if names!="ppfx_nominal":
    #    continue
    #print "Checking this lateral systematics ",names
    nhist = (histo.GetVertErrorBand(names)).GetNHists()
    #print names," has ",nhist," universes"
    hist = histo.GetVertErrorBand(names).GetHist(0)
    hist.GetXaxis().SetTitle("Neutrino Energy GeV")
    hist.GetYaxis().SetTitle("Universe/CV")    
    nbins = hist.GetNbinsX()
    #for i in range(0,100):
       #univfluxval = hist.Interpolate(i*0.1)
       #print "i and value " ,i*0.1,"  ",univfluxval
    hist.Draw("hist")
    hist.SetMaximum(2)
    hist.SetMinimum(0)
    for i in range(0,nhist):
        temphist = (histo.GetVertErrorBand(names)).GetHist(i)
        temphist.Divide(hist)
        temphist.Draw("histeSame")

raw_input()
