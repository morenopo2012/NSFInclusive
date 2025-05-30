import ROOT
from PlotUtils import MnvH1D, MnvH2D, MnvVertErrorBand,MnvLatErrorBand
import os,sys,math
def PrintHistograms(hist,x_axis,y_axis,name):
    canvas = ROOT.TCanvas("c","c",750,750)
    canvas.cd()
    hist.GetXaxis().SetTitle(x_axis)
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().CenterTitle()
    hist.GetYaxis().SetTitle(y_axis)
    hist.Draw("colz")
    raw_input()
    canvas.Print(name+".pdf")
def Get1DAverageChi2(statchi2_hist,NDF=1.0):
    #y axis is the stat universe
    #x axis is the iteration number.
    xbins = statchi2_hist.GetNbinsX()
    ybins = statchi2_hist.GetNbinsY()
    projX_ = statchi2_hist.ProjectionX()
    projX = projX_.Clone()
    projX.Reset()
    for i in range(0,xbins):
        temp_xval = []
        for j in range(0,ybins):
            val = statchi2_hist.GetBinContent(i,j)
            temp_xval.append(val)
        #okay now we have an ensemble of chi2 from stat universe for a given iteraction
        rms = GetRMS(temp_xval)
        avg = sum(temp_xval)/len(temp_xval)
        projX.SetBinContent(i,avg/NDF) #scaling by NDF
        projX.SetBinError(i,rms/NDF)
        projX.GetYaxis().SetTitle("Average #chi ^{2}")
        if NDF>1.0:
            projX.GetYaxis().SetTitle("Average #chi ^{2}/NDF")
        projX.GetYaxis().CenterTitle()
    return projX
def GetRMS(_list):
    mean_val = sum(_list)/len(_list)
    diff_sq = 0.0
    for i in range(0,len(_list)):
        diff_sq += math.pow(mean_val - _list[i],2)*1.0/(float(len(_list)))
    return math.sqrt(diff_sq)
def GetReasonableChi2(old_statchi2_hist):
    statchi2_hist = old_statchi2_hist.Clone("new_statchi2")
    valx=[]
    valy=[]
    max_val=[]
    diff = 1E9
    threshold = 1E3
    while diff>threshold:
        x_index = -99 #iteration number for crazy chi2
        y_index = -99 #stat universe number for crazy chi2
        xbins = statchi2_hist.GetNbinsX() #iterations
        ybins = statchi2_hist.GetNbinsY() #stat universes
        temp_max_val = statchi2_hist.GetMaximum()
        val = []
        for i in range(0,xbins+1):
            for j in range(0,ybins+1):
                temp_val = statchi2_hist.GetBinContent(i,j)
                if temp_val==temp_max_val:
                    x_index=i
                    y_index=j
                    if y_index not in valy:
                        if x_index<10:#no need to worry what happens at higher iteration
                            valx.append(x_index)
                            valy.append(y_index)
                            max_val.append(temp_val)
                else:
                    val.append(temp_val)
        #remove the maximum value calculate the average...
        if len(val)!=0:
            diff = abs(sum(val)/len(val)-temp_max_val)
            #print diff
        else:
            diff = 0.0
        if diff>threshold:
            statchi2_hist.SetBinContent(x_index,y_index,0.0)
    for i in range(0,len(max_val)):
        print " Iteration no.: ",valx[i], " Stat Universe No.: ",valy[i]," Maximum Chi2 :",max_val[i]
    return statchi2_hist
#########################################     
_infile = ROOT.TFile(str(sys.argv[1]),"read")
data = _infile.Get("Input_Hists/h_data")
data_truth = _infile.Get("Input_Hists/h_data_truth")
mc = _infile.Get("Input_Hists/h_mc_reco")
mc_truth = _infile.Get("Input_Hists/h_mc_truth")
migration = _infile.Get("Input_Hists/h_migration_matrix")
#always use modelData_trueData for unfolding test.
stat_chi2 = _infile.Get("Chi2_Iteration_Dists/h_chi2_modelData_trueData_iter_stat")
iter_chi2 = _infile.Get("Chi2_Iteration_Dists/h_chi2_modelData_trueData_iter_chi2")
iter_chi2.GetXaxis().SetRangeUser(1,20)
iter_chi2.GetYaxis().SetRangeUser(0,1000)
c_iter = ROOT.TCanvas("c","c",750,750)
c_iter.cd()
iter_chi2.Draw("colz")
raw_input()
#PrintHistograms(data,"Pz_{#mu} (GeV)","pT_{#mu} (GeV)","data")
#PrintHistograms(data_truth,"Pz_{#mu} (GeV)","pT_{#mu} (GeV)","data_truth")
#PrintHistograms(mc,"Pz_{#mu} (GeV)","pT_{#mu} (GeV)","mc_reco")
#PrintHistograms(mc_truth,"Pz_{#mu} (GeV)","pT_{#mu} (GeV)","mc_truth")
#PrintHistograms(migration,"Bin no","Bin no","migration")
#now the chi2 
new_hist = GetReasonableChi2(stat_chi2) 
c1 = ROOT.TCanvas("c1","c1",750,750)
c1.cd()
#new_hist.GetXaxis().SetRangeUser(0,10)
stat_chi2.Draw("colz")
raw_input()
canvas = ROOT.TCanvas("c","c",750,750)
canvas.cd()
Projhist = Get1DAverageChi2(new_hist,109)  #247 for Q2_enu
Projhist.SetMaximum(5.0)
Projhist.SetMinimum(0.0)
Projhist.GetXaxis().SetRangeUser(1,20)
Projhist.Draw("histe")
Projhist.SaveAs("abc.png")
raw_input()
