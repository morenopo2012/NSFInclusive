#import PlotUtils
import ROOT
from ROOT import PlotUtils
import os,sys

def normalizeHist(obj):

    norm_obj = obj.Clone("%s_row_norm"%(obj.GetName()))
    norm_obj.Reset()

    nbin_y = obj.GetNbinsY()+2
    nbin_x = obj.GetNbinsX()+2

    for i in range(0,nbin_x):
        column_integral = 0
        for j in range(0,nbin_y):
            column_integral += obj.GetBinContent(i,j)
        for j in range(0,nbin_y):
            if column_integral==0: norm_obj.SetBinContent(i,j,0)
            else: norm_obj.SetBinContent(i,j,obj.GetBinContent(i,j)/column_integral)


    return norm_obj

def minmaxHist(obj):

    norm_obj = obj.Clone("%s_row_norm"%(obj.GetName()))
    norm_obj.Reset()

    nbin_y = obj.GetNbinsY()+2
    nbin_x = obj.GetNbinsX()+2

    for i in range(0,nbin_x):
        column_max = 0
        for j in range(0,nbin_y):
            if(obj.GetBinContent(i,j)>column_max): column_max=obj.GetBinContent(i,j)
        for j in range(0,nbin_y):
            if column_max==0: norm_obj.SetBinContent(i,j,0)
            else: norm_obj.SetBinContent(i,j,obj.GetBinContent(i,j)/column_max)


    return norm_obj



inputfile = ROOT.TFile(sys.argv[1])
inputfile.ls()
hist = inputfile.Get(sys.argv[2])
factor = int(sys.argv[3])
hist.Rebin2D(factor,factor)
hist_norm = normalizeHist(hist)
hist_minmax = minmaxHist(hist)

can = ROOT.TCanvas("C1","c1",10,10,1000,700)
hist_norm.Draw("COLZ")
#can.print("RowNormalized.png")

can2 =ROOT.TCanvas("C3","c3",10,10,1000,700)
hist_minmax.Draw("COLZ")
#can2.print("MinMax.png")

raw_input("done")
