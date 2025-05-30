import ROOT
#import PlotUtils
from ROOT import PlotUtils

import os,sys

ROOT.gROOT.SetBatch(True)

inputfile = ROOT.TFile(sys.argv[1])
variable = sys.argv[2]
rebinsize = 2

#myhist = inputfile.Get("h_%sdiff"%(variable))
myhist = inputfile.Get("selected_mc_response2d_%s_migration"%(variable))
myhist.Rebin2D(5,1)
myhist.GetXaxis().SetTitle("%s"%(variable))
nbins_x = myhist.GetNbinsX()
nbins_y = myhist.GetNbinsY()

g1 = ROOT.TGraph()
g1_m = ROOT.TGraph()
g2 = ROOT.TGraph()
g2_m = ROOT.TGraph()

g3 = ROOT.TGraph()
g3_m = ROOT.TGraph()


maxrms = 0
mywidth = 0
#No rebinning
for i in range(1,nbins_x+1):
    myproj = myhist.ProjectionY("raw_%d"%(i),i,i)
    mymaxbin = myproj.GetMaximumBin()
    mymaxbincenter = myproj.GetBinCenter(mymaxbin)
    if(maxrms<myproj.GetRMS()): maxrms=myproj.GetRMS()
    mywidth = myproj.GetBinWidth(i)*10
    myfunc = ROOT.TF1("mygaus","gaus",mymaxbincenter-mywidth,mymaxbincenter+mywidth)
    myproj.Fit(myfunc,"QR+")

    g1.SetPoint(i-1,myhist.GetXaxis().GetBinCenter(i),myfunc.GetParameter(2))
    g1_m.SetPoint(i-1,myhist.GetXaxis().GetBinCenter(i),myfunc.GetParameter(1))

#First rebinning
myhist.Rebin2D(1,rebinsize)
for i in range(1,nbins_x+1):
    myproj = myhist.ProjectionY("raw2_%d"%(i),i,i)
    mymaxbin = myproj.GetMaximumBin()
    mymaxbincenter = myproj.GetBinCenter(mymaxbin)
    myfunc = ROOT.TF1("mygaus","gaus",mymaxbincenter-mywidth,mymaxbincenter+mywidth)
    myproj.Fit(myfunc,"QR+")

    g2.SetPoint(i-1,myhist.GetXaxis().GetBinCenter(i),myfunc.GetParameter(2))
    g2_m.SetPoint(i-1,myhist.GetXaxis().GetBinCenter(i),myfunc.GetParameter(1))

#Second rebinning
myhist.Rebin2D(1,rebinsize)
for i in range(1,nbins_x+1):
    myproj = myhist.ProjectionY("raw3_%d"%(i),i,i)
    mymaxbin = myproj.GetMaximumBin()
    mymaxbincenter = myproj.GetBinCenter(mymaxbin)
    myfunc = ROOT.TF1("mygaus","gaus",mymaxbincenter-mywidth,mymaxbincenter+mywidth)
    myproj.Fit(myfunc,"QR+")

    g3.SetPoint(i-1,myhist.GetXaxis().GetBinCenter(i),myfunc.GetParameter(2))
    g3_m.SetPoint(i-1,myhist.GetXaxis().GetBinCenter(i),myfunc.GetParameter(1))

g2.SetMarkerColor(2)
g2_m.SetMarkerColor(2)
g2.SetMarkerStyle(4)
g2_m.SetMarkerStyle(4)

g3.SetMarkerColor(4)
g3_m.SetMarkerColor(4)
g3.SetMarkerStyle(4)
g3_m.SetMarkerStyle(4)


leg = ROOT.TLegend(0.8,0.8,1,1)
leg.AddEntry(g1,"No Rebin","P")
leg.AddEntry(g2,"1st Rebin","P")
leg.AddEntry(g3,"2nd Rebin","P")


can1 = ROOT.TCanvas("mean","mean",10,10,1000,700)
g1_m.Draw("AP")
g1_m.GetXaxis().SetTitle("%s"%(variable))
g1_m.GetYaxis().SetTitle("Fit Mean")
g2_m.Draw("SAMEP")
g3_m.Draw("SAMEP")
leg.Draw("SAME")


can2 = ROOT.TCanvas("sigma","sigma",10,10,1000,700)
g1.Draw("AP")
g1.GetYaxis().SetTitle("Fit Sigma")
g1.GetXaxis().SetTitle("%s"%(variable))
g1.SetMinimum(0)
g1.SetMaximum(maxrms)
g2.Draw("SAMEP")
g3.Draw("SAMEP")
leg.Draw("SAME")


can1.Print("%s_mean.png"%(variable))
can2.Print("%s_sigma.png"%(variable))



#raw_input("Done")

