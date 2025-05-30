import ROOT
from PlotUtils import MnvH1D,MnvH2D,MnvVertErrorBand
import os,sys
#open the file
if len(sys.argv)<2:
    print "checksystematics <filename>"
    sys.exit(1)
filename  = str(sys.argv[1])
hist_name = str(sys.argv[2])
rootfile = ROOT.TFile("$HISTS/v1/Hists_Energy_rob_v1_t1_z26_Nu_v1_.root","READONLY")
histo = rootfile.Get("Enu_mc").Clone()
vertnames=histo.GetVertErrorBandNames()
            
for names in vertnames:
    #print "Checking this vertical systematics ",names
    nhist = (histo.GetVertErrorBand(names)).GetNHists()
    #if names=="FluxReweight":
    #if names=="Flux":
    if names=="Flux":
        print "It has statuniverse"
    print names," has ",nhist," universes"
    #hist.SetMaximum(2)
    #hist.SetMinimum(0)
    for i in range(0,nhist):
        temphist = (histo.GetVertErrorBand(names)).GetHist(i)
        #temphist.Divide(hist)
        #temphist.Draw(histeSame)
        #integral = temphist.Integral()
        #integral = histo.Integral()
        integral = histo.GetVertErrorBand(names).Integral()
        if integral==0:
            print names,"******* is empty universe (vert)*************"
            break
        #print names, integral 
        #entries = temphist.GetEntries()
        #if entries==0:
        #    print names,"**************** has 0 entries (vert)*******************"
        #    break
print "*****Now printing non empty Systematics*************"
print "*********If name appeared in emtpy systematics******PROBLEM***"           
for names in vertnames:
    #print "Checking this vertical systematics ",names
    nhist = (histo.GetVertErrorBand(names)).GetNHists()
    #print names," has ",nhist," universes"
    #for i in range(0,nhist):
    temphist = (histo.GetVertErrorBand(names)).GetHist(0)
    #integral = temphist.Integral()
    integral = histo.Integral()
    if integral!=0:
        print names," is not  empty(vertlist) ",integral
