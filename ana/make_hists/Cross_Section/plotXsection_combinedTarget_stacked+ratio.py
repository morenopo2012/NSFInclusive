import ROOT
import os,sys
from ROOT import PlotUtils
from ROOT import gStyle, gPad
from ROOT import TLegend, TPad, TColor, TCanvas, TLine

ROOT.gROOT.SetBatch(True)

dirpwd = sys.argv[1]
targetID = sys.argv[2] 
targetZ = sys.argv[3]
plist = sys.argv[4]
scale = sys.argv[5]

mat = None
trueZ = None

if targetZ == "26":
    trueZ = "Iron"
    mat = "Fe"

if targetZ == "82":
    trueZ = "Lead"
    mat = "Pb"

if targetZ == "06":
    trueZ = "Carbon"
    mat = "C"

if targetZ == "99":
    trueZ = "Tracker"
    mat = "CH"

if targetZ == "99":
    #infile= ROOT.TFile(str(dirpwd)+"/CrossSection_t%s_z%02s_%s.root"%(targetID, targetZ, plist))
    infile= ROOT.TFile(str(dirpwd)+"/XSec_thesis_nuclearflux_minervame1L_t99_z99.root")
else:
    infile= ROOT.TFile(str(dirpwd)+"/CrossSection_Daisy_t%s_z%02s_%s.root"%(targetID, targetZ, plist))
#Create a TCanvas on which to draw plots and split it into 2 panels
overall = TCanvas("Data/MC plot", "plot", 800, 800)
top = TPad("Overlay", "Overlay", 0, 0.15, 1, 1)
bottom = TPad("Ratio", "Ratio", 0, 0, 1, 0.15 + 0.116)
#Thou shalt Draw() new TPads lest they be blank!
top.Draw()
bottom.Draw()

mnv = PlotUtils.MnvPlotter()

mcPOT = infile.Get("MCPOT").GetVal()
dataPOT = infile.Get("DataPOT").GetVal()

mcScale = dataPOT/mcPOT
if scale == "1":
    mcScale = 1

#vars = ["Enu","x", "pTmu1D"] #, "ThetamuDeg"]
#vars = ["Enu","x"]
vars = ["Eavailable","q3"]


#steps = ['unfolded','unfolded_effCorrected', 'crossSection', 'crossSection_total']
#steps = ['total_unfolded_effCorrected','crossSection', 'crossSection_total']
steps = ['crossSection', 'crossSection_total']


for step in steps:

    for var in vars:

        data_hist = infile.Get("%s_data_%s"%(step, var)) #crossSection_data_Eavailable

        if step == "unfolded":
            mc_hist = infile.Get("total_efficiency_numerator_%s"%(var))
        elif step == "total_unfolded_effCorrected":
            mc_hist = infile.Get("total_simEventRate_%s"%(var))
            mc_hist_QE = infile.Get("total_simEventRate_QE_%s"%(var))
            mc_hist_RES = infile.Get("total_simEventRate_RES_%s"%(var))
            mc_hist_DIS = infile.Get("total_simEventRate_DIS_%s"%(var))
            mc_hist_Other = infile.Get("total_simEventRate_Other_%s"%(var))
            mc_hist_2p2h = infile.Get("total_simEventRate_2p2h_%s"%(var))
        else:
            mc_hist = infile.Get("simulatedCrossSection_minervame1L")
            #mc_hist = infile.Get("simEventRate_%s_mc_%s"%(step, var))
            #mc_hist_QE = infile.Get("simEventRate_QE_%s_mc_%s"%(step, var))
            #mc_hist_RES = infile.Get("simEventRate_RES_%s_mc_%s"%(step, var))
            #mc_hist_DIS = infile.Get("simEventRate_DIS_%s_mc_%s"%(step, var))
            #mc_hist_Other = infile.Get("simEventRate_Other_%s_mc_%s"%(step, var))
            #mc_hist_2p2h = infile.Get("simEventRate_2p2h_%s_mc_%s"%(step, var))
    
        #mc_hist.PopVertErrorBand("NeutronInelasticExclusives")
        #data_hist.PopVertErrorBand("NeutronInelasticExclusives")

        top.cd()

        mc_hist.GetXaxis().SetLabelColor(ROOT.kWhite)


        if var == "Enu":
            mc_hist.GetXaxis().SetTitle("Antineutrino Energy (GeV)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dE_{#bar{#nu}} (10^{-39} cm^{2}/GeV/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3}/GeV")

        if var == "x":
            #mc_hist.GetXaxis().SetTitle("Bjorken x")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dx_{#bar{#nu}} (10^{-39} cm^{2}/x/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")

        if var == "pTmu1D":
            #mc_hist.GetXaxis().SetTitle("Muon p_{T} (GeV/c)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dp_{T} (10^{-39} cm^{2}/(GeV/c)/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")
       
        if var == "pZmu1D":
            mc_hist.GetXaxis().SetTitle("Muon p_{Z} (GeV/c)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/dp_{||}} (10^{-39} cm^{2}/(GeV/c)/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")

        if var == "ThetamuDeg":
            mc_hist.GetXaxis().SetTitle("Muon Angle (Deg)")
            if step == "crossSection":
                mc_hist.GetYaxis().SetTitle("d#sigma/d#theta_{#mu}} (10^{-39} cm^{2}/Deg/nucleon)")
                gStyle.SetTitleSize(0.05,"y")
            elif step == "crossSection_total":
                mc_hist.GetYaxis().SetTitle("#sigma (10^{-38} cm^{2}/CH)")
                gStyle.SetTitleSize(0.05,"y")
            else:
                mc_hist.GetYaxis().SetTitle("Events #times 10^{3} (norm.)")

        if step == "total_unfolded_effCorrected":
            mc_hist.Scale(1/1E3)
            mc_hist_QE.Scale(1/1E3)
            mc_hist_RES.Scale(1/1E3)
            mc_hist_DIS.Scale(1/1E3)
            mc_hist_Other.Scale(1/1E3)
            mc_hist_2p2h.Scale(1/1E3)
            data_hist.Scale(1/1E3)

        if step == "crossSection":
            mc_hist.Scale(mcScale*1E39)
            #mc_hist_QE.Scale(mcScale*1E39)
            #mc_hist_RES.Scale(mcScale*1E39)
            #mc_hist_DIS.Scale(mcScale*1E39)
            #mc_hist_Other.Scale(mcScale*1E39)
            #mc_hist_2p2h.Scale(mcScale*1E39)
            #data_hist.Scale(1E39)
        
        if step == "crossSection_total":
            mc_hist.Scale(1E38)
            mc_hist_QE.Scale(1E38)
            mc_hist_RES.Scale(1E38)
            mc_hist_DIS.Scale(1E38)
            mc_hist_Other.Scale(1E38)
            mc_hist_2p2h.Scale(1E38)
            data_hist.Scale(1E38)

        #print("DIS: " + str(mc_hist_DIS.Integral()/mc_hist.Integral()))
        #print("RES: " + str(mc_hist_RES.Integral()/mc_hist.Integral()))
        mc_hist.GetXaxis().CenterTitle()
        mc_hist.GetYaxis().CenterTitle()
        #print(data_hist.GetBinContent(data_hist.GetNbinsX())) # last bin
        #print(data_hist.GetBinContent(data_hist.GetNbinsX()+1)) # overflow
        #mc_hist.GetXaxis().SetRangeUser(66,173)

        #mnv.ApplyStyle(1)
        if step == "unfolded":
            mc_hist.Scale(mcScale)
        elif step == "total_unfolded_effCorrected":
            mc_hist.Scale(mcScale)
            mc_hist_QE.Scale(mcScale)
            mc_hist_RES.Scale(mcScale)
            mc_hist_DIS.Scale(mcScale)
            mc_hist_Other.Scale(mcScale)
            mc_hist_2p2h.Scale(mcScale)
        else:
            mc_hist.Scale(1)
            #mc_hist_QE.Scale(1)
            #mc_hist_RES.Scale(1)
            #mc_hist_DIS.Scale(1)
            #mc_hist_Other.Scale(1)
            #mc_hist_2p2h.Scale(1)

        #if var == "Enu":
        #    mc_hist.GetXaxis().SetRangeUser(2, 20)
        if step == "crossSection_total":
            mc_hist.Scale(1, "width")
            mc_hist_QE.Scale(1, "width")
            mc_hist_RES.Scale(1, "width")
            mc_hist_DIS.Scale(1, "width")
            mc_hist_Other.Scale(1, "width")
            mc_hist_2p2h.Scale(1, "width")
            data_hist.Scale(1, "width")
        
        if step == "total_unfolded_effCorrected":
            mc_hist.Scale(1, "width")
            mc_hist_QE.Scale(1, "width")
            mc_hist_RES.Scale(1, "width")
            mc_hist_DIS.Scale(1, "width")
            mc_hist_Other.Scale(1, "width")
            mc_hist_2p2h.Scale(1, "width")
            data_hist.Scale(1, "width")
        #data_hist_stat =  data_hist.GetCVHistoWithStatError() # stat error
        #data_hist_total = data_hist.GetCVHistoWithError() # total error
        #data_hist_sys = data_hist.GetCVHistoWithError(False) # sys error (bool is include stat)
        mc_hist_stat = mc_hist.GetCVHistoWithStatError() 
        
        # MC
        mc_hist.SetLineWidth(3)
        mc_hist.SetLineColor(2)
        # MC error
        mc_hist_stat.SetFillColor(ROOT.kRed-10)
        mc_hist_stat.SetFillStyle(1001)
        mc_hist_stat.SetMarkerStyle(0)
        mc_hist_stat.SetLineWidth(3)
        mc_hist_stat.SetLineColor(2)

        mc_hist.SetLineColor(ROOT.kRed)
        mc_hist.SetLineWidth(2)
        if step == "crossSection_total":
            mc_hist.GetYaxis().SetRangeUser(0, data_hist.GetMaximum()*1.3)
        elif step == "total_unfolded_effCorrected":
            mc_hist.SetMaximum(data_hist.GetMaximum()*1.6)
            if var == "x":
                mc_hist.SetMaximum(data_hist.GetMaximum()*2.0)
                if targetZ == "82":
                    mc_hist.SetMaximum(data_hist.GetMaximum()*2.2)    
        else:
            #mc_hist.SetMaximum(data_hist.GetMaximum()*1.3)
            if var == "x":
                mc_hist.SetMaximum(data_hist.GetMaximum()*1.45)
                if targetZ == "82":
                    mc_hist.SetMaximum(data_hist.GetMaximum()*1.55)  

        # Int channels
        #mc_hist_QE.SetLineWidth(3)
        #mc_hist_QE.SetLineColor(TColor.GetColor("#9A8CE0"))
        #mc_hist_QE.SetFillColor(TColor.GetColor("#9A8CE0"))
        #mc_hist_RES.SetLineColor(TColor.GetColor("#58B778"))
        #mc_hist_RES.SetFillColor(TColor.GetColor("#58B778"))
        #mc_hist_DIS.SetLineColor(TColor.GetColor("#DA8694"))
        #mc_hist_DIS.SetFillColor(TColor.GetColor("#DA8694"))
        #mc_hist_Other.SetLineColor(17)
        #mc_hist_Other.SetFillColor(17)
        #mc_hist_2p2h.SetLineColor(TColor.GetColor("#DDCC77"))
        #mc_hist_2p2h.SetFillColor(TColor.GetColor("#DDCC77"))

        #gStyle.SetErrorX()
        mc_hist.Draw("HIST")
        mc_hist.SetMinimum(0.001)
        #gStyle.SetErrorX(0.5)
        #mc_hist_stat.Draw("E2 SAME")
        #gStyle.SetErrorX(0)
        #mc_hist.Draw("HIST SAME")

        #stack = ROOT.THStack("stack","stack")
        #stack.Add(mc_hist_Other)
        #stack.Add(mc_hist_2p2h)
        #stack.Add(mc_hist_QE)
        #stack.Add(mc_hist_DIS)
        #stack.Add(mc_hist_RES)
        #stack.Draw("HIST SAME")
    
        #gStyle.SetErrorX(0)

        #mc_hist_QE.Draw("HIST SAME")
        #mc_hist_RES.Draw("HIST SAME")
        #mc_hist_DIS.Draw("HIST SAME")
        #mc_hist_Other.Draw("HIST SAME")
        #mc_hist_2p2h.Draw("HIST SAME")
        #data_hist_stat.Draw("SAME HIST p E1 X0")
        #data_hist_stat.Draw("E2 SAME")
        #data_hist_total.Draw("HIST E1 p SAME X0")
        #data_hist.SetMinimum(0.001)

        gPad.RedrawAxis()
        gPad.Update()
    


        ''' void MnvPlotter::DrawDataMCWithErrorBand(
        const MnvH1D* dataHist,
        const MnvH1D* mcHist,
        const Double_t mcScale /*= 1.0*/,
        const std::string& legPos /*= "L"*/,
        const bool useHistTitles /*=false*/,
        const MnvH1D* bkgdHist /*= NULL*/,
        const MnvH1D* dataBkgdHist /*= NULL*/,
        const bool covAreaNormalize/*= false, Area Normalize considerations for covariance matrix*/,
        const bool statPlusSys /* = false */,
        const bool isSmooth /* = false */)
        '''

        #chi2 = mnv.Chi2DataMC(data_hist, mc_hist,mcScale,False,True )
        #ndf = mc_hist.GetNbinsX()
        #print(chi2)
        #if step == "unfolded":
        #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        #elif step == "total_unfolded_effCorrected":
        #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        #elif step == "crossSection_total":
        #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.7, 0.82)
        #else:
        #mnv.AddPOTNormBox(dataPOT,mcPOT, 0.3, 0.82)
        #if var == "x":
        #    mnv.AddPlotLabel("POT normalised", 0.35, 0.75, 0.033, 12, 42)
        #    mnv.AddPlotLabel("Data POT "+ "{:.2E}".format(dataPOT), 0.35, 0.71, 0.033, 12, 42)
        #else:
        #    mnv.AddPlotLabel("POT normalised", 0.55, 0.75, 0.033, 12, 42)
        #    mnv.AddPlotLabel("Data POT "+ "{:.2E}".format(dataPOT), 0.55, 0.71, 0.033, 12, 42)
        if var == "x":
            mnv.AddPlotLabel("POT normalised", 0.65, 0.81, 0.033, 12, 42)
            mnv.AddPlotLabel("Data POT "+ "{:.2E}".format(dataPOT), 0.65, 0.77, 0.033, 12, 42)
        if var == "pTmu1D":
            mnv.AddPlotLabel("POT normalised", 0.65, 0.5, 0.033, 12, 42)
            mnv.AddPlotLabel("Data POT "+ "{:.2E}".format(dataPOT), 0.65, 0.46, 0.033, 12, 42)

        '''
        Chi2DataMC( dataHist, mcHist, ndf, mcScale, useDataErrorMatrix, useOnlyShapeErrors)
        If your data distribution has ONLY STATISTICAL uncertainties, useDataErrorMatrix should be FALSE
        If you are using POT-normalization, then useOnlyShapeErrors should be FALSE
        '''
        if step == "crossSection":
            #chi2 = mnv.Chi2DataMC(data_hist, mc_hist, mcScale, True, False, True)
            # total error from data, statistical error from MC
            ndf = mc_hist.GetNbinsX()
            #print(chi2/ndf)
            # NDF is the number of bins because I am not fitting anything, just comparing rate chi2, refer to L. Fields publications
            # https://arxiv.org/pdf/1305.2234.pdf
            #labelChi2 = "#chi^{2}/ndf = %3.2f/%d = %3.2f"%(chi2, ndf, chi2/ndf)
            #if var == "x":
            #    mnv.AddPlotLabel( labelChi2, .65, .87, 0.033, 12, 42 )
            #else:
            #    mnv.AddPlotLabel( labelChi2, .35, .87, 0.033, 12, 42 )

            
        #if step == 'crossSection':
        #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TL", 0.035, 0.0, True, False)
        
        #elif step == 'crossSection_total':
        #mnv.AddChi2Label(data_hist, mc_hist, mcScale, "TR", 0.035, 0.0, True, False)
        
        if targetZ == "99":
            mnv.AddHistoTitle("Scintillator Cross-Section", 0.045, 62)
        if targetZ == "26":
            mnv.AddHistoTitle("Iron Cross-Section", 0.045, 62)
        if targetZ == "06":
            mnv.AddHistoTitle("Carbon Cross-Section", 0.045, 62)
        if targetZ == "82":
            mnv.AddHistoTitle("Lead Cross-Section", 0.045, 62)

        mc_hist.GetYaxis().SetTitleOffset(0.99)

        if step == "crossSection_total":
            legend = TLegend(0.20,0.45,0.50,0.89)
        elif var == "x":
            legend = TLegend(0.2,0.58,0.5,0.89)
            legend.SetNColumns(1)
        else:
            legend = TLegend(0.55,0.58,0.85,0.89)
        
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.035)
        #legend.AddEntry(mc_hist_stat, " Simulation", "fl")
        #legend.AddEntry(data_hist_stat, " Data", "ep")
        #legend.AddEntry(mc_hist_RES, " RES", "f")
        #legend.AddEntry(mc_hist_DIS, " DIS", "f")
        #legend.AddEntry(mc_hist_QE, " QE", "f")
        #legend.AddEntry(mc_hist_2p2h, " 2p2h", "f")
        #legend.AddEntry(mc_hist_Other, " Other", "f")
        legend.SetTextFont(42)
        legend.Draw()
        top.RedrawAxis()

        top.SetLogx(False)
        if var == "x":
            top.SetLogx()

        bottom.cd()
        bottom.SetTopMargin(0)
        bottom.SetBottomMargin(0.4)
        bottom.SetLogx(False)
        if var == "x":
            bottom.SetLogx()
        #gStyle.SetErrorX(1)


        #ratio = data_hist.GetCVHistoWithStatError().Clone() # stat
        #ratio.Divide(ratio,mc_hist.GetCVHistoWithStatError()) # stat
        #ratio.SetLineWidth(1)
        '''
        for bin in range(0, ratio.GetXaxis().GetNbins()+1):
            print(ratio.GetBinContent(bin))
        if var == "x":
            ratio.GetXaxis().SetTitle("Bjorken x")
        if var == "pTmu1D":
            ratio.GetXaxis().SetTitle("Muon p_{T} (GeV/c)")
        ratio.GetYaxis().SetTitle("Data/MC")
        ratio.GetYaxis().CenterTitle()
        ratio.GetYaxis().SetTitleSize(0.12)
        ratio.GetYaxis().SetTitleOffset(0.3)
        ratio.GetXaxis().CenterTitle()
        ratio.GetYaxis().CenterTitle()
        ratio.GetXaxis().SetTitleOffset(1.15)
        ratio.Draw("X0")

        ratio.SetMaximum(1.5)
        ratio.SetMinimum(0.5)
        ratio.GetXaxis().SetTitleSize(0.16)
        ratio.GetXaxis().SetLabelSize(0.15)
        ratio.GetYaxis().SetLabelSize(0.11)
        ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions
        ratio.GetYaxis().SetTitleSize(0.16)
        ratio.GetXaxis().SetLabelColor(ROOT.kBlack)
        '''

        # Copied from MINERvA-101-Cross-Section/backgroundStack.py
        # same as in void MnvPlotter::DrawDataMCRatio() in MnvPlotter
        # systematic error centered at y = 1
        #sysError = data_hist.GetTotalError(False, True, False) # False for stat error, True for as frac, False for covAreaNorm
        #for whichBin in range(1, sysError.GetXaxis().GetNbins()+1):
        #    sysError.SetBinError(whichBin, max(sysError.GetBinContent(whichBin), 1e-9))
        #    sysError.SetBinContent(whichBin, 1)

        #sysError.SetFillColor(ROOT.kRed-10)
        #sysError.SetFillStyle(1001)
        #sysError.SetMarkerStyle(0)
        #sysError.Draw("E2 SAME")

        # line at 1
        #line = TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
        #line.SetLineColor(46)
        #line.SetLineWidth(2)
        #line.SetLineStyle(9)
        #line.Draw()
        #ratio.Draw("X0")

        #ratio.Draw("X0 SAME")
        gPad.RedrawAxis()
        gPad.Update()

        if targetZ == "99":
            overall.Print("Stacked_%s_t%s_z%02s_%s_%s.png"%(step, targetID, targetZ, var, plist))
            overall.Print("Stacked_%s_t%s_z%02s_%s_%s.pdf"%(step, targetID, targetZ, var, plist))
            overall.Print("Stacked_%s_t%s_z%02s_%s_%s.C"%(step, targetID, targetZ, var, plist))
        else:
            overall.Print("Stacked_%s_Daisy_t%s_z%02s_%s_%s.png"%(step, targetID, targetZ, var, plist))
            overall.Print("Stacked_%s_Daisy_t%s_z%02s_%s_%s.pdf"%(step, targetID, targetZ, var, plist))
            overall.Print("Stacked_%s_Daisy_t%s_z%02s_%s_%s.C"%(step, targetID, targetZ, var, plist))



data_hist.SetDirectory(0)
mc_hist.SetDirectory(0)

print("DONE %s %s %02s"%(plist, targetID, targetZ))
